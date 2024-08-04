#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InverseMap {
    next: Vec<usize>,
}

impl InverseMap {
    pub fn new(size: usize) -> Self {
        let mut next = vec![0; size];
        for i in 0..size {
            next[i] = i;
        }
        Self { next }
    }
    /// merges v0 and v1, must not already be merged.
    pub fn merge(&mut self, v0: usize, v1: usize) {
        //assert!(!self.is_merged(v0, v1));
        let v0_n = self.next[v0];
        let v1_n = self.next[v1];

        self.next[v0] = v1_n;
        self.next[v1] = v0_n;
    }
    pub fn is_merged(&self, v0: usize, v1: usize) -> bool {
        if v0 == v1 {
            return true;
        }
        let mut curr = self.next[v0];
        while curr != v0 {
            if curr == v1 {
                return true;
            }
            curr = self.next[curr];
        }
        return false;
    }
    pub fn merged(&self, v0: usize) -> impl Iterator<Item = usize> + '_ {
        let mut curr = self.next[v0];
        std::iter::repeat(())
            .map(move |()| {
                let n = curr;
                curr = self.next[curr];
                n
            })
            .take_while(move |&c| c != v0)
            .chain(std::iter::once(v0))
    }
}

#[test]
fn test_inv_map() {
    let mut inv_map = InverseMap::new(10);
    assert_eq!(inv_map.merged(0).collect::<Vec<_>>(), vec![0]);
    inv_map.merge(0, 1);
    assert_eq!(inv_map.merged(0).collect::<Vec<_>>(), vec![1, 0]);
    inv_map.merge(2, 3);
    assert_eq!(inv_map.merged(2).collect::<Vec<_>>(), vec![3, 2]);
    inv_map.merge(0, 2);
    assert_eq!(inv_map.merged(2).collect::<Vec<_>>(), vec![1, 0, 3, 2]);
}
