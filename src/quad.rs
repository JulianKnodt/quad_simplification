use super::{cross, dot, kmul, length, normalize, sub, F};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    n: [F; 3],
    p: [F; 3],
}

impl Plane {
    pub fn new_from_quad([v0, v1, v2, v3]: [[F; 3]; 4]) -> Self {
        let n = cross(sub(v0, v2), sub(v1, v3));
        Self::new(v0, n)
    }
    pub fn new(p: [F; 3], n: [F; 3]) -> Self {
        let n = normalize(n);
        Self { n, p }
    }
    /// Point to plane distance
    pub fn dist(&self, p: [F; 3]) -> F {
        let v = sub(self.p, p);
        let perp_v = kmul(dot(v, self.n), self.n);
        length(perp_v)
    }
}
