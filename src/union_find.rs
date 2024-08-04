use std::cell::Cell;
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnionFind {
    ptrs: Vec<Cell<usize>>,

    len: usize,
}

impl UnionFind {
    #[inline]
    pub fn new(size: usize) -> Self {
        let ptrs = vec![Cell::new(0); size];
        for i in 0..size {
            ptrs[i].set(i);
        }
        Self { ptrs, len: size }
    }
    pub fn get(&self, mut v: usize) -> usize {
        while self.ptrs[v].get() != v {
            v = self.ptrs[v].get();
        }
        v
    }
    pub fn get_compress(&self, v: usize) -> usize {
        let dst = self.get(v);
        self.ptrs[v].set(dst);
        dst
    }
    pub fn set(&mut self, v: usize, to: usize) {
        let root_to = self.get_compress(to);
        let root_v = self.get_compress(v);
        if root_v != root_to {
            self.ptrs[root_v].set(root_to);
            self.len -= 1;
        }
    }
    /// Checks if a vertex is itself the root of a tree
    pub fn is_root(&self, v: usize) -> bool {
        self.ptrs[v].get() == v
    }
    pub fn compress(&mut self) {
        for i in 0..self.ptrs.len() {
            // compress it to last item always to flatten pointer chains.
            let terminal = self.get(i);
            if terminal != i {
                self.set(i, terminal);
            }
        }
    }
    pub fn inverse_map(&self) -> Vec<Vec<usize>> {
        let n = self.capacity();
        let mut out = vec![vec![]; n];
        for i in 0..n {
            out[self.get(i)].push(i);
        }
        out
    }
    #[inline]
    pub fn capacity(&self) -> usize {
        self.ptrs.len()
    }
    #[inline]
    pub fn curr_len(&self) -> usize {
        self.len
    }
}
