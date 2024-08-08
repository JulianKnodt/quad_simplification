use super::F;

use std::array::from_fn;

#[inline]
pub fn cross([x, y, z]: [F; 3], [a, b, c]: [F; 3]) -> [F; 3] {
    [y * c - z * b, z * a - x * c, x * b - y * a]
}

#[test]
fn test_cross() {
    assert_eq!(cross([0., 0., 1.], [1., 0., 0.]), [0., 1., 0.]);
    assert_eq!(cross([0., 0., 2.], [1., 0., 0.]), [0., 2., 0.]);
    assert_eq!(cross([1., 0., 0.], [0., 0., 1.]), [0., -1., 0.]);
}

#[inline]
pub fn add<const N: usize>(a: [F; N], b: [F; N]) -> [F; N] {
    from_fn(|i| a[i] + b[i])
}

#[inline]
pub fn sub<const N: usize>(a: [F; N], b: [F; N]) -> [F; N] {
    from_fn(|i| a[i] - b[i])
}

/// L-2 norm of a vector
pub fn length<const N: usize>(v: [F; N]) -> F {
    dot(v, v).sqrt()
}

pub fn norm_inf<const N: usize>(v: [F; N]) -> F {
    v.into_iter().map(F::abs).max_by(F::total_cmp).unwrap()
}

pub fn to_3<const N: usize>(v: [F; N]) -> [F; 3] {
    assert!(N >= 3);
    [v[0], v[1], v[2]]
}

pub fn zero_cat([v0, v1, v2]: [F; 3]) -> [F; 4] {
    [v0, v1, v2, 0.]
}

pub fn one_cat([v0, v1, v2]: [F; 3]) -> [F; 4] {
    [v0, v1, v2, 1.]
}

#[inline]
pub fn normalize<const N: usize>(v: [F; N]) -> [F; N] {
    let sum: F = v.iter().map(|v| v * v).sum();
    if sum < 1e-10 {
        return [0.; N];
    }
    let s = sum.sqrt().recip();
    v.map(|v| v * s)
}

#[inline]
pub fn dot<const N: usize>(a: [F; N], b: [F; N]) -> F {
    (0..N).map(|i| a[i] * b[i]).sum()
}

#[inline]
pub fn kmul<const N: usize>(k: F, xyz: [F; N]) -> [F; N] {
    xyz.map(|v| v * k)
}

#[inline]
pub fn minmax<const N: usize>(vs: impl Iterator<Item = [F; N]>) -> [[F; N]; 2] {
    vs.fold([[F::INFINITY; N], [F::NEG_INFINITY; N]], |[min, max], n| {
        use std::array::from_fn;
        [from_fn(|i| min[i].min(n[i])), from_fn(|i| max[i].max(n[i]))]
    })
}

/// Computes the conjugate for inverse rotation of a quaternion.
#[inline]
pub fn conj([x, y, z, w]: [F; 4]) -> [F; 4] {
    [-x, -y, -z, w]
}

fn quat_mul([r1, r2, r3, r0]: [F; 4], [s1, s2, s3, s0]: [F; 4]) -> [F; 4] {
    [
        r0 * s1 + r1 * s0 - r2 * s3 + r3 * s2,
        r0 * s2 + r1 * s3 + r2 * s0 - r3 * s1,
        r0 * s3 - r1 * s2 + r2 * s1 + r3 * s0,
        r0 * s0 - r1 * s1 - r2 * s2 - r3 * s3,
    ]
}

pub fn quat_rot([x, y, z]: [F; 3], quat: [F; 4]) -> [F; 3] {
    let v = [x, y, z, 0.];
    let [a, b, c, _] = quat_mul(quat_mul(quat, v), conj(quat));
    [a, b, c]
}

#[inline]
pub fn tri_area(a: [F; 3], b: [F; 3], c: [F; 3]) -> F {
    length(cross(sub(a, b), sub(a, c))) / 2.
}

/// Computes the quad area given a set of vertices.
#[inline]
pub fn quad_area([q0, q1, q2, q3]: [[F; 3]; 4]) -> F {
    tri_area(q0, q1, q2) + tri_area(q0, q2, q3)
    //0.5 * length(sub(q1, q3)) * length(sub(q0, q2))
}

#[inline]
pub fn poly_area(mut v: impl Iterator<Item = [F; 3]>) -> F {
    let mut acc = 0.;
    let Some(fst) = v.next() else {
        return 0.;
    };
    let Some(mut v0) = v.next() else {
        return 0.;
    };
    for v1 in v {
        acc += tri_area(fst, v0, v1);
        v0 = v1;
    }
    acc
}

#[inline]
pub fn quat_from_to(s: [F; 3], t: [F; 3]) -> [F; 4] {
    let v = cross(t, s);
    normalize([v[0], v[1], v[2], 1. + dot(normalize(s), normalize(t))])
}

#[inline]
pub fn quat_from_axis_angle(axis: [F; 3], angle: F) -> [F; 4] {
    let s = (angle / 2.).sin();
    let [x, y, z] = axis.map(|v| v * s);
    [x, y, z, (angle / 2.).cos()]
}

/// Computes rotation from the standard xyz basis to this basis, where fwd and up are orthogonal
/// and normalized.
pub fn quat_from_standard(fwd: [F; 3], up: [F; 3]) -> [F; 4] {
    let r0 = quat_from_to([1., 0., 0.], fwd);
    let r1 = quat_from_to(quat_rot([0., 1., 0.], r0), up);
    quat_mul(r1, r0)
}

#[test]
fn test_quat() {
    let q = quat_from_to([1., 0., 0.], [0., 1., 0.]);
    let rot = quat_rot([1., 0., 0.], q);
    assert!(length(sub(rot, [0., 1., 0.])) < 1e-3);
}

#[test]
pub fn test_quat_basis() {
    let tgt = normalize([0., 0.5, 0.5]);
    let up = [1., 0., 0.];

    let q = quat_from_standard(tgt, up);

    let r0 = quat_rot([1., 0., 0.], q);
    let r1 = quat_rot([0., 1., 0.], q);
    assert!(length(sub(r0, tgt)) < 1e-4);
    assert!(length(sub(r1, up)) < 1e-4);
}
