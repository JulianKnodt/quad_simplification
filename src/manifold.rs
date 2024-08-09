use crate::inv_map::InverseMap;
use crate::union_find::UnionFind;

use std::cmp::minmax;

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum FaceKind {
    Degenerate,
    Tri([usize; 3]),
    Quad([usize; 4]),
    Polygon(Vec<usize>),
}

impl FaceKind {
    pub fn from_slice(s: &[usize]) -> Self {
        match s {
            &[] | &[_] | &[_, _] => Self::Degenerate,
            &[a, b, _] | &[a, _, b] | &[_, a, b] if a == b => Self::Degenerate,

            &[a, b, c] => {
                let mut t = [a, b, c];
                t.sort_unstable();
                Self::Tri(t)
            }
            &[a, b, c, d]
            | &[d, a, b, c]
            | &[c, d, a, b]
            | &[a, c, d, b]
            | &[a, c, b, d]
            | &[d, a, c, b]
                if a == b =>
            {
                Self::from_slice(&[a, c, d])
            }
            &[a, b, c, d] => {
                let mut q = [a, b, c, d];
                q.sort_unstable();
                Self::Quad(q)
            }
            v => Self::Polygon(v.to_vec()),
        }
    }
    pub fn has_edge(&self, [e0, e1]: [usize; 2]) -> bool {
        let s = self.as_slice();
        (0..s.len()).any(|i| {
            (s[i] == e0 && s[(i + 1) % s.len()] == e1) || (s[i] == e1 && s[(i + 1) % s.len()] == e0)
        })
    }
    pub fn remap(&mut self, map: impl Fn(usize) -> usize) {
        match self {
            FaceKind::Degenerate => return,
            FaceKind::Tri(vs) => {
                let mut new = vs.map(map);
                new.sort_unstable();
                let [a, b, c] = new;
                *self = if a == b || b == c {
                    FaceKind::Degenerate
                } else {
                    FaceKind::Tri(new)
                }
            }
            FaceKind::Quad(vs) => {
                let mut new = vs.map(map);
                new.sort_unstable();
                *self = match new {
                    // only two possibilities because it's sorted
                    [a, b, c, d] | [d, a, b, c] if (a == b && b == c) || (a == b && c == d) => {
                        FaceKind::Degenerate
                    }
                    [a, b, c, d] | [d, a, b, c] | [c, d, a, b] if a == b => {
                        FaceKind::Tri([a, c, d])
                    }
                    new => FaceKind::Quad(new),
                }
            }
            _ => todo!(),
        }
    }
    pub fn is_degenerate(&self) -> bool {
        matches!(self, FaceKind::Degenerate)
    }
    /// Returns the set of vertices in this face as a single slice.
    pub fn as_slice(&self) -> &[usize] {
        match self {
            FaceKind::Degenerate => &[],
            FaceKind::Tri(vs) => vs.as_slice(),
            FaceKind::Quad(vs) => vs.as_slice(),
            FaceKind::Polygon(vs) => &vs,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum EdgeKind {
    #[default]
    Empty,
    Boundary(usize),
    Manifold([usize; 2]),
    NonManifold(Vec<usize>),
}

impl EdgeKind {
    pub fn len(&self) -> usize {
        match self {
            EdgeKind::Empty => 0,
            EdgeKind::Boundary(_) => 1,
            EdgeKind::Manifold(_) => 2,
            EdgeKind::NonManifold(ref fs) => fs.len(),
        }
    }
    pub fn insert(&mut self, face: usize) {
        match self {
            &mut EdgeKind::Empty => {
                *self = EdgeKind::Boundary(face);
            }
            &mut EdgeKind::Boundary(f0) => {
                assert_ne!(f0, face);
                *self = EdgeKind::Manifold([f0, face]);
            }
            &mut EdgeKind::Manifold([f0, f1]) => {
                assert_ne!(f0, face);
                assert_ne!(f1, face);
                *self = EdgeKind::NonManifold(vec![f0, f1, face]);
            }
            EdgeKind::NonManifold(ref mut fs) => {
                assert!(!fs.contains(&face));
                fs.push(face);
            }
        }
    }
    pub fn as_slice(&self) -> &[usize] {
        match self {
            EdgeKind::Empty => &[],
            EdgeKind::Boundary(v) => core::slice::from_ref(v),
            EdgeKind::Manifold(vs) => vs,
            EdgeKind::NonManifold(vs) => vs.as_slice(),
        }
    }
    pub fn is_empty(&self) -> bool {
        self.as_slice().is_empty()
    }

    /// Deletes equivalent faces between two EdgeKinds, and moves them into `self`.
    pub fn drain(&mut self, o: &mut Self) {
        let o: Self = std::mem::replace(o, EdgeKind::Empty);
        let n = match (&self, o) {
            (EdgeKind::Boundary(a), EdgeKind::Boundary(b)) if *a == b => EdgeKind::Empty,
            (EdgeKind::Boundary(_), EdgeKind::Boundary(_)) => todo!(),

            (EdgeKind::Manifold([a, b]), EdgeKind::Boundary(c))
            | (EdgeKind::Manifold([b, a]), EdgeKind::Boundary(c))
                if c == *b =>
            {
                EdgeKind::Boundary(*a)
            }

            (EdgeKind::Manifold(_), EdgeKind::Boundary(_)) => todo!(),

            (EdgeKind::Manifold([a, b]), EdgeKind::Manifold([c, d]))
            | (EdgeKind::Manifold([a, b]), EdgeKind::Manifold([d, c]))
            | (EdgeKind::Manifold([b, a]), EdgeKind::Manifold([c, d]))
            | (EdgeKind::Manifold([b, a]), EdgeKind::Manifold([d, c]))
                if *b == c =>
            {
                EdgeKind::Manifold([*a, d])
            }

            (_, o) => todo!("{self:?} {o:?}"),
        };
        *self = n;
    }
    pub fn pairwise(&self) -> impl Iterator<Item = [usize; 2]> + '_ {
        let v = self.as_slice();
        (0..v.len()).flat_map(move |i| ((i + 1)..v.len()).map(move |j| [v[i], v[j]]))
    }
}

/// A mesh representation which is suitable for collapsing vertices.
/// It can associate data with each vertex.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CollapsibleManifold<T> {
    vertices: UnionFind,

    pub edges: Vec<Vec<usize>>,

    inv_map: InverseMap,

    pub data: Vec<T>,
}

impl<T: Copy> CollapsibleManifold<T> {
    #[inline]
    pub fn new(size: usize) -> Self
    where
        T: Default,
    {
        Self {
            vertices: UnionFind::new(size),
            edges: vec![],
            inv_map: InverseMap::new(size),

            data: vec![T::default(); size],
        }
    }
    pub fn new_with(size: usize, f: impl Fn(usize) -> T) -> Self {
        let mut data = Vec::with_capacity(size);
        for i in 0..size {
            data.push(f(i));
        }
        Self {
            vertices: UnionFind::new(size),
            edges: vec![vec![]; size],
            inv_map: InverseMap::new(size),

            data,
        }
    }
    pub fn num_vertices(&self) -> usize {
        self.vertices.curr_len()
    }
    pub fn vertices(&self) -> impl Iterator<Item = (usize, &T)> + '_ {
        (0..self.vertices.capacity())
            .filter(|&vi| !self.deleted(vi))
            .map(|vi| (vi, &self.data[vi]))
    }
    /// Computes all vertices merged into other vertices
    pub fn source_vertices(&self) -> Vec<Vec<usize>> {
        self.vertices.inverse_map()
    }
    pub fn deleted(&self, vi: usize) -> bool {
        !self.vertices.is_root(vi)
    }

    /// Adds an edge. For faces, should call `add_face`.
    pub fn add_edge(&mut self, v0: usize, v1: usize) {
        if v0 == v1 {
            return;
        }
        let [v0, v1] = minmax(v0, v1);
        self.edges[v0].push(v1);
        self.edges[v1].push(v0);
    }

    pub fn add_face(&mut self, face: &[usize]) {
        let n = face.len();
        for i in 0..n {
            self.add_edge(face[i], face[(i + 1) % n]);
        }
    }

    pub fn vertex_adj(&self, v: usize) -> impl Iterator<Item = usize> + '_ {
        self.edges[v]
            .iter()
            .map(|&dst| self.vertices.get_compress(dst))
    }

    pub fn merge(&mut self, v0: usize, v1: usize, merge: impl Fn(T, &mut T)) {
        assert_ne!(v0, v1);
        let [src, dst] = minmax(v0, v1);
        assert!(!self.deleted(src));
        assert!(!self.deleted(dst));

        self.vertices.set(src, dst);

        self.inv_map.merge(src, dst);

        merge(self.data[src], &mut self.data[dst]);
        self.data[src] = self.data[dst];

        // correctly set adjacent faces for each edge below
        let mut src_e = std::mem::take(&mut self.edges[src]);
        self.edges[dst].append(&mut src_e);
        self.edges[dst].retain(|&e1| self.vertices.get_compress(e1) != dst);
        self.edges[dst].sort_unstable_by_key(|e1| self.vertices.get_compress(*e1));
        self.edges[dst].dedup_by_key(|e1| self.vertices.get(*e1));
    }
    pub fn merged_vertices(&self, v0: usize) -> impl Iterator<Item = usize> + '_ {
        self.inv_map.merged(v0)
    }

    pub fn get(&self, v: usize) -> &T {
        &self.data[self.vertices.get(v)]
    }
    pub fn set(&mut self, v: usize, t: T) {
        self.data[self.vertices.get(v)] = t;
    }

    /// All edges in this manifold mesh with v0-v1 in sorted order.
    pub fn ord_edges(&self) -> impl Iterator<Item = [usize; 2]> + '_ {
        self.edges.iter().enumerate().flat_map(move |(src, dsts)| {
            dsts.iter()
                .map(move |&dst| minmax(src, self.vertices.get(dst)))
        })
    }
}
