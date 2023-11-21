use crate::assign::AssignedPoint;
use crate::circuit::ecc_chip::{EccChipBaseOps, EccChipScalarOps};
use crate::circuit::integer_chip::IntegerChipOps;
use crate::context::Context;
use crate::context::GeneralScalarEccContext;
use crate::tests::{random_bls12_381_fr, run_circuit_on_bn256};
use crate::utils::field_to_bn;
use ark_std::{end_timer, start_timer};
use halo2_proofs::pairing::bls12_381::{G1Affine, G1, Fr as Fb};
use halo2_proofs::pairing::bn256::Fr;
use halo2_proofs::poly::EvaluationDomain;
use std::cell::RefCell;
use std::rc::Rc;

#[test]
fn test_bls12_381_polyeval_chip_over_bn256_fr() {
    let ctx = Rc::new(RefCell::new(Context::new()));
    let mut ctx = GeneralScalarEccContext::<G1Affine, Fr>::new(ctx);

    let values = vec![1, 2, 3, 4];

    let values:Vec<Fb> = values
        .clone()
        .into_iter()
        .map(|x| Fb::from(x))
        .collect();

    let x = random_bls12_381_fr();

    // for _ in 0..50 {
    //     let a = random_bls12_381_fr();
    //     let b = random_bls12_381_fr();
    //     let p = G1::generator() * a;
    //     acc = acc + p * b;
    //     points.push(p);
    //     scalars.push(b);
    // }

    let n = values.len(); 

    let logn = (n as f32).log2() as u32;

    let domain = EvaluationDomain::<Fb>::new(2, logn);

    let omega = domain.get_omega();    

    let mut invs = Vec::with_capacity(4096);

    let mut evaldomain = Vec::with_capacity(4096);

    let mut acc = Fb::from(0);

    let mut omega_i = Fb::from(1);

    let mut x_n = Fb::from(1);
    
    for i in 0..n {
        evaldomain.push(omega_i);

        let inv = (x - omega_i).invert().unwrap();
        invs.push(inv);

        acc += (values[i]) * omega_i * inv;

        omega_i = omega_i * omega;

        x_n = x_n * x;
        //print!("outside acc value:{:?}", acc);
    }

    acc = (x_n - Fb::one()) * Fb::from(n as u64).invert().unwrap() * acc;


    let timer = start_timer!(|| "setup");

    let values = values
        .into_iter()
        .map(|x| ctx.scalar_integer_ctx.assign_w(&field_to_bn(&x)))
        .collect::<Vec<_>>();

    let evaldomain = evaldomain
        .into_iter()
        .map(|x| ctx.scalar_integer_ctx.assign_w(&field_to_bn(&x)))
        .collect::<Vec<_>>(); 



    let x = ctx.scalar_integer_ctx.assign_w(&field_to_bn(&x));

    let res = ctx.eval(&values, &x, &evaldomain);

    let res_expect = ctx.scalar_integer_ctx.assign_w(&field_to_bn(&acc));

    ctx.scalar_integer_ctx.is_int_equal(&res, &res_expect);

    end_timer!(timer);

    run_circuit_on_bn256(ctx.into(), 22);
}