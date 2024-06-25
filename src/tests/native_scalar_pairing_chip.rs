use crate::assign::{AssignedCondition, AssignedG2Affine, self};
use crate::circuit::base_chip::BaseChipOps;
use crate::circuit::ecc_chip::{EccChipBaseOps, EccChipScalarOps, EccBaseIntegerChipWrapper};
use crate::circuit::fq12::{Fq12ChipOps, Fq2ChipOps};
use crate::circuit::pairing_chip::PairingChipOps;
use crate::context::IntegerContext;
use crate::context::{Context, NativeScalarEccContext};
use crate::tests::run_circuit_on_bn256;
use crate::utils::field_to_bn;
use halo2_proofs::arithmetic::{CurveAffine, Field};
use halo2_proofs::pairing::bn256::pairing;
use halo2_proofs::pairing::bn256::{Fr, G1Affine, G2Affine, G2, Gt};
use halo2_proofs::pairing::group::cofactor::CofactorCurveAffine;
use rand::rngs::OsRng;
use std::cell::RefCell;
use std::ops::Mul;
use std::rc::Rc;

use super::bench_circuit_on_bn256;

fn build_bn256_pairing_chip_over_bn256_fr_circuit() -> NativeScalarEccContext<G1Affine> {
    // {
    //     let ctx = Rc::new(RefCell::new(Context::new()));
    //     let ctx = IntegerContext::<halo2_proofs::pairing::bn256::Fq, Fr>::new(ctx);
    //     let mut ctx = NativeScalarEccContext::<G1Affine>(ctx);

    //     let a = G1Affine::random(&mut OsRng);
    //     let b = G2Affine::from(G2::random(&mut OsRng));

    //     let ab = pairing(&a, &b);

    //     let bx = ctx.fq2_assign_constant((
    //         b.coordinates().unwrap().x().c0,
    //         b.coordinates().unwrap().x().c1,
    //     ));
    //     let by = ctx.fq2_assign_constant((
    //         b.coordinates().unwrap().y().c0,
    //         b.coordinates().unwrap().y().c1,
    //     ));
    //     let b = AssignedG2Affine::new(
    //         bx,
    //         by,
    //         AssignedCondition(ctx.0.ctx.borrow_mut().assign_constant(Fr::zero())),
    //     );

    //     let ab0 = ctx.fq12_assign_constant((
    //         (
    //             (ab.0.c0.c0.c0, ab.0.c0.c0.c1),
    //             (ab.0.c0.c1.c0, ab.0.c0.c1.c1),
    //             (ab.0.c0.c2.c0, ab.0.c0.c2.c1),
    //         ),
    //         (
    //             (ab.0.c1.c0.c0, ab.0.c1.c0.c1),
    //             (ab.0.c1.c1.c0, ab.0.c1.c1.c1),
    //             (ab.0.c1.c2.c0, ab.0.c1.c2.c1),
    //         ),
    //     ));

    //     let a = ctx.assign_point(&a.to_curve());

    //     let ab1 = ctx.pairing(&[(&a, &b)]);

    //     ctx.fq12_assert_eq(&ab0, &ab1);

    //     ctx
    //     // println!("run circuit");

    //     // run_circuit_on_bn256(ctx.into(), 20);

    // }
    
    //our public input will be
    //CT and Tk = {C_i} {D_i} C'    {w_i} L' {R'_i}  R'  =>  T
    //public input circuit : instance  hash =>  {C_i} {D_i} C'    {w_i}  L' {R'_i}  R'   (2*i+1)*64 bytes  32*i bytes  (i+2) * 128bytes = 4i+2+4i+8+i = 9i+10  i=100 910*32bytes 
    {
        //suppose we have G1: C_i, G2: L', G1:D_i,  G2:R'_i  scalar：w_i
        //then we can compute
        // C_w_i = ctx.ecc_mul(C_i, w_i)  D_w_i = ctx.ecc_mul(D_i, w_i) 
        // T_i = ctx.pairing ((C_w_i,L'],(D_w_i, R'i))
        // ctx.fq12_assert_eq  fq12_mul(T , \prod T_i) = ctx.pairing(C', R')
        // 目标 输出T  证明T是通过G1: C_i, G2: L', G1:D_i,  G2:R'_i  scalar：w_i计算得出的
        let i = 2;
        let ctx = Rc::new(RefCell::new(Context::new()));
        let ctx = IntegerContext::<halo2_proofs::pairing::bn256::Fq, Fr>::new(ctx);
        let mut ctx = NativeScalarEccContext::<G1Affine>(ctx);

        let w = vec![Fr::random(&mut OsRng); i];
        let C = vec![G1Affine::random(&mut OsRng); i];
        let D = vec![G1Affine::random(&mut OsRng); i];
        let L = G2Affine::from(G2::random(&mut OsRng));
        let R = vec![G2Affine::from(G2::random(&mut OsRng)); i];

        // let C_W_1 = G1Affine::from(C_1.mul(w_1));
        let C_W = w.iter().zip(C.iter()).map(|(wi,ci)| G1Affine::from(ci.mul(wi))).collect::<Vec<_>>();

        let D_W = w.iter().zip(D.iter()).map(|(wi,di)| G1Affine::from(di.mul(wi))).collect::<Vec<_>>();
        
        let mut T = pairing(&C_W[0], &L) + pairing(&D_W[0], &R[0]);

        for index in 1..i {
            T += pairing(&C_W[index], &L) + pairing(&D_W[index], &R[index]);
        }

        let real_T = ctx.fq12_assign_constant((
            (
                (T.0.c0.c0.c0, T.0.c0.c0.c1),
                (T.0.c0.c1.c0, T.0.c0.c1.c1),
                (T.0.c0.c2.c0, T.0.c0.c2.c1),
            ),
            (
                (T.0.c1.c0.c0, T.0.c1.c0.c1),
                (T.0.c1.c1.c0, T.0.c1.c1.c1),
                (T.0.c1.c2.c0, T.0.c1.c2.c1),
            ),
        ));


        //in circuit

        let assigned_C = C.iter().map(|c| ctx.assign_point(&c.to_curve())).collect::<Vec<_>>();

        let assigned_D = D.iter().map(|d| ctx.assign_point(&d.to_curve())).collect::<Vec<_>>();
        
        let assigned_w = w.iter().map(|w|ctx.0.ctx.borrow_mut().assign_constant(*w)).collect::<Vec<_>>();

        let assigned_C_w = assigned_C.iter().zip(assigned_w.iter()).map(|(c, w)|ctx.ecc_mul(&c, *w)).collect::<Vec<_>>();

        let assigned_D_w = assigned_D.iter().zip(assigned_w.iter()).map(|(d, w)|ctx.ecc_mul(&d, *w)).collect::<Vec<_>>();

        let Lx = ctx.fq2_assign_constant((
            L.coordinates().unwrap().x().c0,
            L.coordinates().unwrap().x().c1,
        ));
        let Ly = ctx.fq2_assign_constant((
            L.coordinates().unwrap().y().c0,
            L.coordinates().unwrap().y().c1,
        ));
        let assigned_L = AssignedG2Affine::new(
            Lx,
            Ly,
            AssignedCondition(ctx.0.ctx.borrow_mut().assign_constant(Fr::zero())),
        );
        let assigned_R = R.iter().map(|r|
            AssignedG2Affine::new(
                ctx.fq2_assign_constant((
                    r.coordinates().unwrap().x().c0,
                    r.coordinates().unwrap().x().c1,
                )),
                ctx.fq2_assign_constant((
                    r.coordinates().unwrap().y().c0,
                    r.coordinates().unwrap().y().c1,
                )),
                AssignedCondition(ctx.0.ctx.borrow_mut().assign_constant(Fr::zero())),
            )
        ).collect::<Vec<_>>();


        let mut terms = vec![];

        for index in 0..i {
            terms.push((&assigned_C_w[index], &assigned_L));
            terms.push((&assigned_D_w[index], &assigned_R[index]));
        }


        let recover_T = ctx.pairing(&terms);

        ctx.fq12_assert_eq(&real_T, &recover_T);
        ctx
        
    }

    // {
    //     let ctx = Rc::new(RefCell::new(Context::new()));
    //     let ctx = IntegerContext::<halo2_proofs::pairing::bn256::Fq, Fr>::new(ctx);
    //     let mut ctx = NativeScalarEccContext::<G1Affine>(ctx);

    //     let a = G1Affine::random(&mut OsRng);
    //     let b = G2Affine::from(G2::random(&mut OsRng));

    //     let bx = ctx.fq2_assign_constant((
    //         b.coordinates().unwrap().x().c0,
    //         b.coordinates().unwrap().x().c1,
    //     ));
    //     let by = ctx.fq2_assign_constant((
    //         b.coordinates().unwrap().y().c0,
    //         b.coordinates().unwrap().y().c1,
    //     ));
    //     let b = AssignedG2Affine::new(
    //         bx,
    //         by,
    //         AssignedCondition(ctx.0.ctx.borrow_mut().assign_constant(Fr::zero())),
    //     );

    //     let neg_a = ctx.assign_point(&-a.to_curve());
    //     let a = ctx.assign_point(&a.to_curve());

    //     ctx.check_pairing(&[(&a, &b), (&neg_a, &b)]);

    //     ctx
    // }
}

#[test]
fn test_bn256_pairing_chip_over_bn256_fr() {
    let ctx = build_bn256_pairing_chip_over_bn256_fr_circuit();
    run_circuit_on_bn256(ctx.into(), 22);
}

#[test]
fn bench_bn256_pairing_chip_over_bn256_fr() {
    let ctx = build_bn256_pairing_chip_over_bn256_fr_circuit();
    bench_circuit_on_bn256(ctx.into(), 22);
}
