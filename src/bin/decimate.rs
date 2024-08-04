#![feature(cmp_minmax)]
#![feature(array_windows)]
use clap::Parser;
use pars3d::obj;

use priority_queue::PriorityQueue;

use quad_collapse::manifold::{CollapsibleManifold, EdgeKind, FaceKind};
use quad_collapse::{add, cross, kmul, normalize, poly_area, sub, F};

use ordered_float::NotNan;
use std::collections::HashMap;

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


    let mut pq: PriorityQueue<EdgeCollapse, NotNan<F>> =
        PriorityQueue::with_capacity(m.num_vertices());
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
