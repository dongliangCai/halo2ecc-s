use crate::assign::AssignedPoint;
use crate::circuit::base_chip::BaseChipOps;
use crate::circuit::ecc_chip::{EccChipBaseOps, EccChipScalarOps};
use crate::context::IntegerContext;
use crate::context::{Context, NativeScalarEccContext};
use crate::tests::{random_fr, run_circuit_on_bn256};
use ark_std::{end_timer, start_timer};
use halo2_proofs::pairing::bn256::{Fr, G1Affine, G1};
use halo2_proofs::pairing::group::Group;
use std::cell::RefCell;
use std::rc::Rc;

#[test]
fn test_native_ecc_chip() {
    let ctx = Rc::new(RefCell::new(Context::new()));
    let ctx = IntegerContext::<halo2_proofs::pairing::bn256::Fq, Fr>::new(ctx);
    let mut ctx = NativeScalarEccContext(ctx, 0);

    let mut points = vec![];
    let mut scalars = vec![];
    let mut acc = G1::identity();

    for _ in 0..100 {
        let a = random_fr();
        let b = random_fr();
        let p = G1::generator() * a;
        acc = acc + p * b;
        points.push(p);
        scalars.push(b);
    }

    let timer = start_timer!(|| "setup");

    let points = points
        .iter()
        .map(|x| ctx.assign_point(x))
        .collect::<Vec<_>>();
    let scalars = scalars
        .into_iter()
        .map(|x| ctx.0.ctx.borrow_mut().assign(x))
        .collect::<Vec<_>>();
    let res: AssignedPoint<_, _> = ctx.msm(&points, &scalars);
    let res_expect: AssignedPoint<G1Affine, Fr> = ctx.assign_point(&acc);
    ctx.ecc_assert_equal(&res, &res_expect);

    end_timer!(timer);

    run_circuit_on_bn256(ctx.into(), 22);
}

#[test]
fn test_native_ecc_unsafe_chip() {
    let ctx = Rc::new(RefCell::new(Context::new()));
    let ctx = IntegerContext::<halo2_proofs::pairing::bn256::Fq, Fr>::new(ctx);
    let mut ctx = NativeScalarEccContext::<G1Affine>(ctx, 0);

    let mut points = vec![];
    let mut scalars = vec![];
    let mut acc = G1::identity();

    for _ in 0..100 {
        let a = random_fr();
        let b = random_fr();
        let p = G1::generator() * a;
        acc = acc + p * b;
        points.push(p);
        scalars.push(b);
    }

    let timer = start_timer!(|| "setup");

    let points = points
        .iter()
        .map(|x| ctx.assign_nonzero_point(x))
        .collect::<Vec<_>>();
    let scalars = scalars
        .into_iter()
        .map(|x| ctx.0.ctx.borrow_mut().assign(x))
        .collect::<Vec<_>>();

    let res = ctx.msm_batch_on_group_unsafe(
        &points,
        &scalars,
        G1::generator() * random_fr(),
        G1::generator() * random_fr(),
    );
    let res_expect = ctx.assign_nonzero_point(&acc);
    ctx.ecc_assert_equal_unsafe(&res, &res_expect);

    end_timer!(timer);

    run_circuit_on_bn256(ctx.into(), 22);
}
