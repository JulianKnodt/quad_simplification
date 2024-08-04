/// An index into a vector of arrays.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash, Default)]
struct VecOfArraysIdx {
    num_vertices: u8,
    idx: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VecOfArrays<T> {
    data: Vec<Vec<T>>,

    num_elems: usize,
}

/// reinterpret a vector of elements as a vector of arrays.
unsafe fn reinterpret<T, const N: usize>(v: Vec<T>) -> Vec<[T; N]> {
    let (ptr, len, cap) = v.into_raw_parts();
    assert_eq!(len % N, 0);
    assert_eq!(cap % N, 0);
    unsafe { Vec::from_raw_parts(ptr as *mut [T; N], len / N, cap / N) }
}

impl<T> VecOfArrays<T> {
    pub fn new() -> Self {
        Self {
            data: vec![],
            num_elems: 0,
        }
    }
    pub fn push<const N: usize>(&mut self, arr: [T; N]) {
        if self.data.len() <= N {
            self.data.resize_with(N + 1, Vec::new);
        }
        self.num_elems += 1;

        let mut v = unsafe { reinterpret(std::mem::take(&mut self.data[N])) };
        v.push(arr);
        self.data[N] = v.into_flattened();
    }
    pub fn push_slice(&mut self, slice: &[T])
    where
        T: Copy,
    {
        let n = slice.len();
        if self.data.len() <= n {
            self.data.resize_with(n + 1, Vec::new);
        }
        self.num_elems += 1;

        self.data[n].extend_from_slice(slice);
        // have to be careful here to ensure that capacity is valid for [T; N]
        // explicitly round it up so that it is divisible by n
        let cap = self.data[n].capacity();
        self.data[n].shrink_to(cap + (n - (cap % n)));
    }
    pub fn get<const N: usize>(&self) -> &[[T; N]] {
        if self.data.len() <= N {
            return &[];
        }
        unsafe { self.data[N].as_chunks_unchecked::<N>() }
    }
    pub fn max_n(&self) -> usize {
        self.data.len()
    }
    pub fn len_n<const N: usize>(&self) -> usize {
        self.data[N].len() / N
    }
    pub fn len(&self) -> usize {
        self.num_elems
    }
    /// For getting chunks at run-time
    pub fn get_chunks(&self, n: usize) -> impl Iterator<Item = &[T]> + '_ {
        self.data
            .get(n)
            .into_iter()
            .flat_map(move |d| d.chunks_exact(n))
    }

    pub fn get_idx(&self, mut idx: usize) -> &[T] {
        let mut n = 0;
        while idx >= self.data[n].len() {
            idx -= self.data[n].len();
            n += 1;
        }
        self.get_chunks(n).nth(idx).unwrap()
    }
    pub fn update_mut<const N: usize>(&mut self, f: impl Fn(&mut Vec<[T; N]>)) {
        if self.data.len() <= N {
            return;
        }
        let mut v = unsafe { reinterpret(std::mem::take(&mut self.data[N])) };
        f(&mut v);
        self.data[N] = v.into_flattened();
    }
}

#[test]
fn test_basic() {
    let mut v = VecOfArrays::new();
    v.push([0, 1, 2]);
    assert_eq!(&[[0, 1, 2]], v.get::<3>());
    v.update_mut::<3>(|x| assert_eq!(x.as_slice(), &[[0, 1, 2]]));
}
