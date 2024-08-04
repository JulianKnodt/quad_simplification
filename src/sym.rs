use std::ops::{Add, Mul, Neg, Sub};

use super::F;

pub const fn sum_up_to<const N: usize>() -> usize {
    N * (N + 1) / 2
}

fn sqr(x: F) -> F {
    x * x
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SymMatrix<const N: usize>
where
    [(); sum_up_to::<N>()]:,
{
    data: [F; sum_up_to::<N>()],
}

pub type SymMatrix4 = SymMatrix<4>;
pub type SymMatrix3 = SymMatrix<3>;

macro_rules! define_det {
  ($N: expr, $p: vis $name: ident, $name2: ident,
    $x0: ident,
    [ $( $op: tt $x: ident, $y: ident -> [$($a: ident, $b: ident),+]),+],
    $next_det: ident, $next_det2: ident) => {

    #[inline]
    $p fn $name <$(const $x: usize, const $y: usize,)+>(&self) -> F {
      let mut acc = 0.;
      $(
        let m = self.v_const::<$x0, $y>();
        acc $op m * self.$next_det::< $($a, $b,)+ >();
      )+
      acc
    }

    #[inline]
    $p fn $name2 <const COL: usize, $(const $x: usize, const $y: usize,)+>(&self, col: [F; $N]) -> F {
      let mut acc = 0.;
      $(
        let m = self.v_with_col_const::<COL, $x0, $y>(col);
        acc $op m * self.$next_det2::<COL, $($a, $b,)+ >(col);
      )+
      acc
    }
  }
}
impl<const N: usize> SymMatrix<N>
where
    [(); sum_up_to::<N>()]:,
{
    #[inline]
    pub fn new(data: [F; sum_up_to::<N>()]) -> Self {
        Self { data }
    }
    pub fn zero() -> Self {
        Self { data: [0.; _] }
    }
    pub fn is_zero(&self) -> bool {
        self.data.iter().all(|&v| v == 0.)
    }
}

impl SymMatrix4 {
    #[rustfmt::skip]
    pub const SYM_IDX: [[usize; 4];4] = [
      [0,1,2,3],
      [1,4,5,6],
      [2,5,7,8],
      [3,6,8,9],
    ];
    #[inline]
    pub const fn v(&self, x: usize, y: usize) -> F {
        self.data[Self::SYM_IDX[x][y]]
    }

    /// Constant version of access
    #[inline]
    pub const fn v_const<const X: usize, const Y: usize>(&self) -> F {
        self.data[Self::SYM_IDX[X][Y]]
    }

    #[inline]
    pub const fn v_with_col_const<const COL: usize, const X: usize, const Y: usize>(
        &self,
        col: [F; 4],
    ) -> F {
        if COL == X {
            col[Y]
        } else {
            self.v_const::<X, Y>()
        }
    }

    pub fn col(&self, j: usize) -> [F; 4] {
        std::array::from_fn(|i| self.v(j, i))
    }

    #[inline]
    pub fn qr(&self) -> ([[F; 4]; 4], [[F; 4]; 4]) {
        use super::{dot, kmul, length, sub};
        // [col][row]
        let mut v = [[0.; 4]; 4];
        let mut q = [[0.; 4]; 4];
        let mut r = [[0.; 4]; 4];

        for i in 0..4 {
            v[i] = self.col(i);
        }

        for i in 0..4 {
            r[i][i] = length(v[i]);
            if r[i][i] != 0. {
                q[i] = kmul(1. / r[i][i], v[i]);
            }
            for j in i + 1..4 {
                r[i][j] = dot(q[i], v[j]);
                v[j] = sub(v[j], kmul(r[i][j], q[i]));
            }
        }

        (q, r)
    }

    fn eigen_2x2<const X0: usize, const Y0: usize, const X1: usize, const Y1: usize>(
        &self,
    ) -> ([F; 2], [[F; 2]; 2]) {
        let a = self.v_const::<X0, Y0>();
        let b = self.v_const::<X1, Y1>();
        let c = self.v_const::<X0, Y1>();

        if c.abs() < 1e-8 {
            return ([a, b], [[1., 0.], [0., 1.]]);
        }

        let delta = 4. * sqr(c) + sqr(a - b);
        let delta = delta.sqrt();
        debug_assert!(delta.is_finite());

        let [e0, e1] = [delta, -delta].map(|d| (a + b - d) / 2.);

        use super::normalize;
        let vs = [e0, e1].map(|e| normalize([(e - b) / c, 1.]));
        /*
        let prods = vs.map(|v| [a * v[0] + c * v[1], c * v[0] + b * v[1]]);
        assert_eq!(
            super::kmul(e0, vs[0]),
            prods[0],
            "{e0} {:?} [[{a}, {c}], [{c}, {b}]]",
            vs[0]
        );
        assert_eq!(
            super::kmul(e1, vs[1]),
            prods[1],
            "{e1} {:?} [[{a}, {c}], [{c}, {b}]]",
            vs[1]
        );
        */
        ([e0, e1], vs)
    }

    #[inline]
    pub fn eigen_3x3(&self) -> ([F; 3], [[F; 3]; 3]) {
        let p1 =
            sqr(self.v_const::<0, 1>()) + sqr(self.v_const::<0, 2>()) + sqr(self.v_const::<1, 2>());
        let diag = [
            self.v_const::<0, 0>(),
            self.v_const::<1, 1>(),
            self.v_const::<2, 2>(),
        ];
        const ZERO_EPS: F = 1e-8;
        // diagonal
        if p1.abs() < ZERO_EPS {
            return (diag, [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]);
        }
        let tr = diag.into_iter().sum::<F>();
        let q = tr / 3.;
        let p2 = sqr(self.v_const::<0, 0>() - q)
            + sqr(self.v_const::<1, 1>() - q)
            + sqr(self.v_const::<2, 2>() - q)
            + 2. * p1;
        let p = (p2 / 6.).sqrt();

        let b = (*self - SymMatrix4::ident() * q) * (1. / p);
        let r = b.det_qem() / 2.;

        const PI: F = std::f64::consts::PI as F;
        let phi = if r <= -1. {
            PI
        } else if r >= 1. {
            0.
        } else {
            r.acos()
        } / 3.;
        // TODO may need to fix these eigenvalues.

        let eig0 = q + 2. * p * phi.cos();
        let eig2 = q + 2. * p * (phi + 2. * PI / 3.).cos();
        let eig1 = 3. * q - eig0 - eig2;
        let eigs = [eig0, eig1, eig2];

        // https://hal.science/hal-01501221/document
        let a = self.v_const::<0, 0>();
        let b = self.v_const::<1, 1>();
        let c = self.v_const::<2, 2>();
        let d = self.v_const::<0, 1>();
        let e = self.v_const::<1, 2>();
        let f = self.v_const::<0, 2>();
        // TODO ensure that this paper is correct, and also verify eigenvalues are correct
        // FIXME handle degenerate cases as well.

        assert!(d != 0. || e != 0. || f != 0.);
        let (e, v) = match f.abs() < ZERO_EPS {
            true if d.abs() < ZERO_EPS => {
                let ([e0, e1], [v0, v1]) = self.eigen_2x2::<1, 1, 2, 2>();
                let vecs = [[0., v0[0], v0[1]], [0., v1[0], v1[1]], [1., 0., 0.]];
                ([e0, e1, a], vecs)
            }
            true if e.abs() < ZERO_EPS => {
                let ([e0, e1], [v0, v1]) = self.eigen_2x2::<0, 0, 1, 1>();
                let vecs = [[v0[0], v0[1], 0.], [v1[0], v1[1], 0.], [0., 0., 1.]];
                ([e0, e1, c], vecs)
            }
            true => self.to_3x3().power_eigen(),
            false => {
                let denoms = eigs.map(|eig| f * (b - eig) - d * e);
                if denoms.iter().any(|&d| d.abs() < ZERO_EPS) {
                    self.to_3x3().power_eigen()
                } else {
                    let numers = eigs.map(|eig| d * (c - eig) - e * f);
                    let ms: [F; 3] = std::array::from_fn(|i| numers[i] / denoms[i]);
                    let vecs = std::array::from_fn(|i| {
                        super::normalize([(eigs[i] - c - e * ms[i]) / f, ms[i], 1.])
                    });
                    (eigs, vecs)
                }
                /*
                let vecs = std::array::from_fn(|i| {
                    let eig = eigs[i];
                    if eig.abs() < ZERO_EPS {
                        return [0., 0., 0.];
                    }

                    let denom = f * (b - eig) - d * e;
                    let numer = d * (c - eig) - e * f;
                    assert!(
                        denom.abs() > ZERO_EPS,
                        "{denom} {eig} {numer} {eigs:?} {:?}",
                        self.flat::<3>()
                    );
                    let m = numer / denom;
                    super::normalize([(eig - c - e * m) / f, m, 1.])
                });
                (eigs, vecs)
                */
            }
        };

        // check that it is an actual eigen vector
        /*
        for i in 0..3 {
            use super::{kmul, length, sub};
            let m = self.to_3x3().vec_mul(v[i]);
            let o = kmul(e[i], v[i]);
            assert!(length(sub(m, o)) < 1e-4, "{m:?} {o:?}");
        }
        */
        (e, v)
    }

    #[inline]
    pub fn outer([x, y, z, w]: [F; 4]) -> Self {
        Self::new([
            x * x,
            x * y,
            x * z,
            x * w,
            y * y,
            y * z,
            y * w,
            z * z,
            z * w,
            w * w,
        ])
    }

    #[rustfmt::skip]
    pub fn ident() -> Self {
        Self::new([
          1.,0.,0.,0.,
             1.,0.,0.,
                1.,0.,
                   1.
        ])
    }
    #[inline]
    pub fn fill(&mut self, v: F) {
        self.data.fill(v)
    }

    // TODO decide if it's better to implement LU or do this explicit version?
    fn det2x2<const X0: usize, const Y0: usize, const X1: usize, const Y1: usize>(&self) -> F {
        self.v_const::<X0, Y0>() * self.v_const::<X1, Y1>()
            - self.v_const::<X0, Y1>() * self.v_const::<X1, Y0>()
    }
    fn det2x2_with_col<
        const COL: usize,
        const X0: usize,
        const Y0: usize,
        const X1: usize,
        const Y1: usize,
    >(
        &self,
        col: [F; 4],
    ) -> F {
        self.v_with_col_const::<COL, X0, Y0>(col) * self.v_with_col_const::<COL, X1, Y1>(col)
            - self.v_with_col_const::<COL, X0, Y1>(col) * self.v_with_col_const::<COL, X1, Y0>(col)
    }

    define_det!(4, det3x3, det3x3_with_col, X0, [
      +=X0, Y0 -> [X1, Y1, X2, Y2],
      -=X1, Y1 -> [X1, Y0, X2, Y2],
      +=X2, Y2 -> [X1, Y0, X2, Y1]
    ], det2x2, det2x2_with_col);

    define_det!(4, det4x4, det4x4_with_col, X0, [
      +=X0, Y0 -> [X1, Y1, X2, Y2, X3, Y3],
      -=X1, Y1 -> [X1, Y0, X2, Y2, X3, Y3],
      +=X2, Y2 -> [X1, Y0, X2, Y1, X3, Y3],
      -=X3, Y3 -> [X1, Y0, X2, Y1, X3, Y2]
    ], det3x3, det3x3_with_col);

    #[inline]
    pub fn det(&self) -> F {
        self.det4x4::<0, 0, 1, 1, 2, 2, 3, 3>()
    }

    #[inline]
    pub fn det_with_col<const COL: usize>(&self, col: [F; 4]) -> F {
        self.det4x4_with_col::<COL, 0, 0, 1, 1, 2, 2, 3, 3>(col)
    }

    #[inline]
    pub fn det_qem(&self) -> F {
        self.det3x3::<0, 0, 1, 1, 2, 2>()
    }

    #[inline]
    pub fn det_qem_with_col<const COL: usize>(&self, col: [F; 4]) -> F {
        self.det3x3_with_col::<COL, 0, 0, 1, 1, 2, 2>(col)
    }

    #[inline]
    pub fn vec_mul(&self, v: [F; 4]) -> [F; 4] {
        std::array::from_fn(|i| (0..4).map(|j| self.v(i, j) * v[i]).sum())
    }

    pub fn vec_mul3(&self, v: [F; 3]) -> [F; 3] {
        std::array::from_fn(|i| (0..3).map(|j| self.v(i, j) * v[i]).sum())
    }

    /// Matmul two symmetric matrices
    pub fn matmul(&self, o: &Self) -> Self {
        let e = |i, j| (0..4).map(|k| self.v(i, k) * o.v(k, j)).sum();
        Self::new([
            e(0, 0),
            e(1, 0),
            e(2, 0),
            e(3, 0),
            e(1, 1),
            e(2, 1),
            e(3, 1),
            e(2, 2),
            e(3, 2),
            e(3, 3),
        ])
    }

    pub fn flat<const N: usize>(&self) -> [[F; N]; N] {
        std::array::from_fn(|i| std::array::from_fn(|j| self.v(i, j)))
    }

    pub fn to_3x3(&self) -> SymMatrix3 {
        let v = |x, y| self.v(x, y);
        SymMatrix3::new([v(0, 0), v(0, 1), v(0, 2), v(1, 1), v(1, 2), v(2, 2)])
    }

    pub fn square(&self) -> Self {
        self.matmul(self)
    }
}

impl SymMatrix3 {
    #[rustfmt::skip]
    pub fn ident() -> Self {
        Self::new([
          1.,0.,0.,
             1.,0.,
                1.,
        ])
    }

    #[inline]
    pub fn outer([x, y, z]: [F; 3]) -> Self {
        Self::new([x * x, x * y, x * z, y * y, y * z, z * z])
    }
    #[rustfmt::skip]
    pub const SYM_IDX: [[usize;3];3] = [
      [0,1,2],
      [1,3,4],
      [2,4,5],
    ];

    #[inline]
    pub fn v(&self, x: usize, y: usize) -> F {
        self.data[Self::SYM_IDX[x][y]]
    }

    /// Constant version of access
    #[inline]
    pub const fn v_const<const X: usize, const Y: usize>(&self) -> F {
        self.data[Self::SYM_IDX[X][Y]]
    }

    #[inline]
    pub const fn v_with_col_const<const COL: usize, const X: usize, const Y: usize>(
        &self,
        col: [F; 3],
    ) -> F {
        if COL == X {
            col[Y]
        } else {
            self.v_const::<X, Y>()
        }
    }

    pub fn col(&self, x: usize) -> [F; 3] {
        std::array::from_fn(|i| self.v(i, x))
    }
    pub fn nz_col(&self) -> Option<[F; 3]> {
        for i in 0..3 {
            let c = self.col(i);
            if super::length(c) > 1e-3 {
                return Some(c);
            }
        }
        None
    }

    /// Matmul two symmetric matrices
    pub fn matmul(&self, o: &Self) -> Self {
        let e = |i, j| (0..3).map(|k| self.v(i, k) * o.v(k, j)).sum();
        Self::new([e(0, 0), e(1, 0), e(2, 0), e(1, 1), e(2, 1), e(2, 2)])
    }

    pub fn square(&self) -> Self {
        self.matmul(self)
    }

    pub fn vec_mul(&self, v: [F; 3]) -> [F; 3] {
        std::array::from_fn(|i| (0..3).map(|j| self.v(j, i) * v[j]).sum())
    }

    pub fn flat(&self) -> [[F; 3]; 3] {
        std::array::from_fn(|i| std::array::from_fn(|j| self.v(i, j)))
    }
    pub fn power_eigen(&self) -> ([F; 3], [[F; 3]; 3]) {
        let mut q0 = [1., 0., 0.];
        let mut q1 = [0., 1., 0.];
        use super::{cross, dot, kmul, length, normalize, sub};
        for _ in 0..16 {
            let n_q0 = normalize(self.vec_mul(q0));
            if length(n_q0) > 1e-4 {
                q0 = n_q0;
            }
            let n_q1 = self.vec_mul(q1);
            let n_q1 = normalize(sub(n_q1, kmul(dot(n_q1, q0), q0)));
            if length(n_q1) > 1e-4 {
                q1 = n_q1;
            }
        }

        let q = [q0, q1, cross(q0, q1)];
        let eigs = q.map(|col| dot(self.vec_mul(col), col));
        (eigs, q)
    }

    fn det2x2<const X0: usize, const Y0: usize, const X1: usize, const Y1: usize>(&self) -> F {
        self.v_const::<X0, Y0>() * self.v_const::<X1, Y1>()
            - self.v_const::<X0, Y1>() * self.v_const::<X1, Y0>()
    }

    fn det2x2_with_col<
        const COL: usize,
        const X0: usize,
        const Y0: usize,
        const X1: usize,
        const Y1: usize,
    >(
        &self,
        col: [F; 3],
    ) -> F {
        self.v_with_col_const::<COL, X0, Y0>(col) * self.v_with_col_const::<COL, X1, Y1>(col)
            - self.v_with_col_const::<COL, X0, Y1>(col) * self.v_with_col_const::<COL, X1, Y0>(col)
    }

    define_det!(3, det3x3, det3x3_with_col, X0, [
      +=X0, Y0 -> [X1, Y1, X2, Y2],
      -=X1, Y1 -> [X1, Y0, X2, Y2],
      +=X2, Y2 -> [X1, Y0, X2, Y1]
    ], det2x2, det2x2_with_col);

    #[inline]
    pub fn det(&self) -> F {
        self.det3x3::<0, 0, 1, 1, 2, 2>()
    }

    pub fn det_with_col<const COL: usize>(&self, col: [F; 3]) -> F {
        self.det3x3_with_col::<COL, 0, 0, 1, 1, 2, 2>(col)
    }
}

#[test]
fn test_power_eigen() {
    let n = SymMatrix3::new([1., 2., 3., 4., 5., 6.]);
    let (e, _) = n.power_eigen();
    let near = |a: F, b: F| (a - b).abs() < 1e-3;
    println!("{:?}", n.flat());
    assert!(near(e[0], 11.344), "{e:?}");
    assert!(near(e[1], -0.515), "{e:?}");
    assert!(near(e[2], 0.170), "{e:?}");
    println!("{e:?}");
}

impl<const N: usize> Add for SymMatrix<N>
where
    [(); sum_up_to::<N>()]:,
{
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(std::array::from_fn(|i| self.data[i] + o.data[i]))
    }
}

impl<const N: usize> Sub for SymMatrix<N>
where
    [(); sum_up_to::<N>()]:,
{
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(std::array::from_fn(|i| self.data[i] - o.data[i]))
    }
}

impl<const N: usize> Mul<F> for SymMatrix<N>
where
    [(); sum_up_to::<N>()]:,
{
    type Output = Self;
    fn mul(self, o: F) -> Self {
        Self::new(std::array::from_fn(|i| self.data[i] * o))
    }
}

impl<const N: usize> Neg for SymMatrix<N>
where
    [(); sum_up_to::<N>()]:,
{
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(self.data.map(Neg::neg))
    }
}

#[test]
fn test_ident_det() {
    let v = SymMatrix4::ident();
    //assert_eq!(v.det2x2::<0,0,1,1>(), 1.);
    //assert_eq!(v.det3x3::<0,0,1,1,2,2>(), 1.);
    assert_eq!(v.det4x4::<0, 0, 1, 1, 2, 2, 3, 3>(), 1.);
}

#[test]
fn test_diag() {
    #[rustfmt::skip]
    let v = SymMatrix4::new([
      2., 0., 0., 0.,
          1., 0., 0.,
              1., 0.,
                  1.
    ]);
    assert_eq!(v.det(), 2.);

    #[rustfmt::skip]
    let v = SymMatrix4::new([
      1., 0., 0., 0.,
          2., 0., 0.,
              1., 0.,
                  1.
    ]);
    assert_eq!(v.det(), 2.);

    #[rustfmt::skip]
    let v = SymMatrix4::new([
      1., 0., 0., 0.,
          1., 0., 0.,
              2., 0.,
                  1.
    ]);
    assert_eq!(v.det(), 2.);

    #[rustfmt::skip]
    let v = SymMatrix4::new([
      1., 0., 0., 0.,
          1., 0., 0.,
              1., 0.,
                  2.
    ]);
    assert_eq!(v.det(), 2.);
}

#[test]
fn test_det() {
    #[rustfmt::skip]
    let v = SymMatrix4::new([
      1., 2., 0., 0.,
          1., 0., 0.,
              1., 0.,
                  1.
    ]);
    assert_eq!(v.det(), -3.);

    #[rustfmt::skip]
    let v = SymMatrix4::new([
      1., 0., 2., 0.,
          1., 0., 0.,
              1., 0.,
                  1.
    ]);
    assert_eq!(v.det(), -3.);

    #[rustfmt::skip]
    let v = SymMatrix4::new([
      1., 0., 0., 2.,
          1., 0., 0.,
              1., 0.,
                  1.
    ]);
    assert_eq!(v.det(), -3.);

    #[rustfmt::skip]
    let v = SymMatrix4::new([
      1., 0., 0., 0.,
          1., 2., 0.,
              1., 0.,
                  1.
    ]);
    assert_eq!(v.det(), -3.);

    #[rustfmt::skip]
    let v = SymMatrix4::new([
      1., 0., 0., 0.,
          1., 0., 2.,
              1., 0.,
                  1.
    ]);
    assert_eq!(v.det(), -3.);

    #[rustfmt::skip]
    let v = SymMatrix4::new([
      1., 0., 0., 0.,
          1., 0., 0.,
              1., 2.,
                  1.
    ]);
    assert_eq!(v.det(), -3.);
}

#[test]
fn test_eig() {
    let sym = SymMatrix4::new([
        8.0,
        0.0,
        0.0,
        1.1920929e-7,
        0.0,
        0.0,
        0.0,
        8.0,
        -2.3841858e-7,
        8.0,
    ]);
    let ([e0, e1, e2], _) = sym.eigen_3x3();
    assert_eq!(8., e0);
    assert_eq!(0., e1);
    assert_eq!(8., e2);

    let ([e0, e1, e2], _) = SymMatrix4::ident().eigen_3x3();
    assert_eq!(e0, 1.);
    assert_eq!(e1, 1.);
    assert_eq!(e2, 1.);
}

#[test]
fn test_matmul() {
    let v = SymMatrix4::outer([4., 3., 2., 1.]);
    assert_eq!(v.data, [16., 12., 8., 4., 9., 6., 3., 4., 2., 1.]);
}
