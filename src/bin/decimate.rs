#![feature(cmp_minmax)]
#![feature(array_windows)]
use clap::Parser;
use pars3d::obj;

use priority_queue::PriorityQueue;

use quad_collapse::manifold::{CollapsibleManifold, EdgeKind};
use quad_collapse::{add, cross, kmul, normalize, poly_area, quad_area, sub, tri_area, F};

use ordered_float::NotNan;

use std::cmp::minmax;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

#[derive(Parser)]
struct Args {
    #[arg(short, long, required = true)]
    input: String,

    #[arg(short, long, required = true)]
    output: String,

    /// Ratio of faces
    #[arg(short, long, required = true)]
    ratio: F,
}

fn main() {
    let args = Args::parse();

    let mut input_obj = obj::parse(&args.input, false, false)
        .unwrap_or_else(|_| panic!("Failed to parse input obj file: {}", args.input));

    assert_eq!(input_obj.objects.len(), 1);
    let mesh = input_obj.objects.pop().unwrap();

    let target_faces = args.ratio * mesh.f.len() as F;
    let target_faces = target_faces as usize;

    let mut edge_adj = HashMap::<[usize; 2], EdgeKind>::new();

    let mut vert_face = vec![vec![]; mesh.v.len()];

    for (fi, f) in mesh.f.iter().enumerate() {
        for &[v0, v1] in f.v.array_windows::<2>() {
            let e = std::cmp::minmax(v0, v1);
            edge_adj.entry(e).or_default().insert(fi);
        }
        let e = std::cmp::minmax(*f.v.first().unwrap(), *f.v.last().unwrap());
        edge_adj.entry(e).or_default().insert(fi);

        for &vi in &f.v {
            vert_face[vi].push(fi);
        }
    }

    let mut m = CollapsibleManifold::<[usize; 4]>::new_with(mesh.f.len(), |fi| {
        assert_eq!(mesh.f[fi].v.len(), 4, "Temporary Limitation");
        match mesh.f[fi].v.as_slice() {
            &[a, b, c, d] => [a, b, c, d],
            _ => todo!(),
        }
    });

    for vf in vert_face {
        m.add_face(&vf);
    }

    let mut pq: PriorityQueue<FaceMerge, NotNan<F>> =
        PriorityQueue::with_capacity(m.num_vertices());

    macro_rules! add_edge {
        ($f0: expr, $f1: expr) => {{
            let [f0, f1] = minmax($f0, $f1);
            assert!(!m.deleted(f0));
            assert!(!m.deleted(f1));

            let q0 = *m.get(f0);
            let q1 = *m.get(f1);

            let Some((new_quad @ [a, b, x, y], [s0, s1])) = merged_quad(q0, q1) else {
                continue;
            };

            let original_total_area = m
                .merged_vertices(f0)
                .chain(m.merged_vertices(f1))
                .map(|fi| poly_area(mesh.f[fi].v.iter().map(|&vi| mesh.v[vi])))
                .sum::<F>();

            let new_quad_area = quad_area(new_quad.map(|q| mesh.v[q]));

            let scaffold_tris = [[b, s0, y], [a, s1, x]];
            let scaffold_area = scaffold_tris
                .into_iter()
                .map(|vis| {
                    let [v0, v1, v2] = vis.map(|vi| mesh.v[vi]);
                    tri_area(v0, v1, v2)
                })
                .sum::<F>();

            let cost = (original_total_area) - (new_quad_area + scaffold_area);
            let cost = cost.abs();
            assert!(cost.is_finite());

            let face_merge = FaceMerge { f0, f1 };
            pq.push(face_merge, NotNan::new(-cost).unwrap());
        }};
    }

    for [f0, f1] in m.ord_edges() {
        add_edge!(f0, f1);
    }

    while let Some((next, _cost)) = pq.pop() {
        if m.num_vertices() <= target_faces {
            break;
        }

        let FaceMerge { f0, f1 } = next;
        if m.deleted(f0) || m.deleted(f1) {
            continue;
        }

        let q0 = *m.get(f0);
        let q1 = *m.get(f1);

        let Some((nq, _shared)) = merged_quad(q0, q1) else {
            continue;
        };

        // also need to add scaffolding triangles here?

        m.merge(f0, f1, |_, _| nq);
        assert!(m.deleted(f0));
        assert!(!m.deleted(f1));
        for adj in m.vertex_adj(f1) {
            add_edge!(adj, f1);
        }
    }

    let out = File::create(args.output).unwrap();
    let mut out = BufWriter::new(out);
    // write out obj file
    for [x, y, z] in &mesh.v {
        writeln!(out, "v {x} {y} {z}").unwrap();
    }

    for (_, xyzw) in m.vertices() {
        let [x, y, z, w] = xyzw.map(|i| i + 1);
        writeln!(out, "f {x} {y} {z} {w}").unwrap();
    }
}

/// Computes the quad merged together, along with the shared edge
pub fn merged_quad(
    [a, b, c, d]: [usize; 4],
    [x, y, z, w]: [usize; 4],
) -> Option<([usize; 4], [usize; 2])> {
    match ([a, b, c, d], [x, y, z, w]) {
        (
            [a, b, s0, s1] | [s1, a, b, s0] | [s0, s1, a, b] | [b, s0, s1, a],
            [x, y, t0, t1] | [t1, x, y, t0] | [t0, t1, x, y] | [y, t0, t1, x],
        ) if minmax(s0, s1) == minmax(t0, t1) => Some(([a, b, x, y], minmax(t0, t1))),

        _ => None,
    }
}

/// A merge between two faces
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
struct FaceMerge {
    f0: usize,
    f1: usize,
}

impl core::hash::Hash for FaceMerge {
    #[inline]
    fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
        self.f0.hash(state);
        self.f1.hash(state);
    }
}
