#![feature(cmp_minmax)]
#![feature(array_windows)]
use clap::Parser;
use pars3d::obj;

use priority_queue::PriorityQueue;

use quad_collapse::manifold::{CollapsibleManifold, EdgeKind};
use quad_collapse::{poly_area, quad_area, tri_area, F};

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

    // map from each edge to a scaffolds on that edge if it exists
    let mut scaffolds: HashMap<[usize; 2], [usize; 3]> = HashMap::new();
    // TODO also need to account for boundary edges to not scaffold for them.

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

            let scaf0 = scaffolds.get(&minmax(b, s0));
            if scaf0 != scaffolds.get(&minmax(s0, x)) {
                continue;
            }

            let scaf1 = scaffolds.get(&minmax(a, s1));
            if scaf1 != scaffolds.get(&minmax(s1, y)) {
                continue;
            }

            let scaf_area = |scaf: Option<_>, vis: [usize; 3]| {
                let [v0, v1, v2] = vis.map(|vi| mesh.v[vi]);
                let area = tri_area(v0, v1, v2);
                if scaf.is_some() {
                    -area
                } else {
                    area
                }
            };

            let scaffold_area = scaf_area(scaf0, [b, s0, x]) + scaf_area(scaf1, [a, s1, y]);

            let original_total_area = m
                .merged_vertices(f0)
                .chain(m.merged_vertices(f1))
                .map(|fi| poly_area(mesh.f[fi].v.iter().map(|&vi| mesh.v[vi])))
                .sum::<F>();

            let new_quad_area = quad_area(new_quad.map(|q| mesh.v[q]));

            // need to add area of new scaffolds, and delete area of removed scaffolds.
            let cost = original_total_area - new_quad_area - scaffold_area;
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

        let Some((nq @ [a, b, x, y], [s0, s1])) = merged_quad(q0, q1) else {
            continue;
        };

        let scaf0 = scaffolds.get(&minmax(b, s0));
        if scaf0 != scaffolds.get(&minmax(s0, x)) {
            continue;
        }
        let scaf0 = scaf0.copied();

        let scaf1 = scaffolds.get(&minmax(a, s1));
        if scaf1 != scaffolds.get(&minmax(s1, y)) {
            continue;
        }
        let scaf1 = scaf1.copied();

        let scaf0_e = [minmax(b, s0), minmax(s0, x), minmax(x, b)];
        for e in scaf0_e {
            match scaf0 {
                None => assert!(scaffolds.insert(e, [b, s0, x]).is_none()),
                Some(_) => assert!(scaffolds.remove(&e).is_some()),
            }
        }

        let scaf1_e = [minmax(a, s1), minmax(s1, y), minmax(y, a)];
        for e in scaf1_e {
            match scaf1 {
                None => assert!(scaffolds.insert(e, [a, s1, y]).is_none()),
                Some(_) => assert!(scaffolds.remove(&e).is_some()),
            }
        }

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

    for (&[e0, e1], vis) in scaffolds.iter() {
        if minmax(e0, e1) != minmax(vis[0], vis[1]) {
            continue;
        }
        let [i, j, k] = vis.map(|vi| vi + 1);
        writeln!(out, "f {i} {j} {k}").unwrap();
    }
    /*
     */
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
        ) if minmax(s0, s1) == minmax(t0, t1) => {
            assert_ne!(a, t0);
            assert_ne!(b, t0);
            assert_ne!(x, t0);
            assert_ne!(y, t0);

            assert_ne!(a, t1);
            assert_ne!(b, t1);
            assert_ne!(x, t1);
            assert_ne!(y, t1);

            Some(([a, b, x, y], [s0, s1]))
        }

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
