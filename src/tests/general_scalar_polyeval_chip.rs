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
use halo2_proofs::pairing::group::ff::PrimeField;
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

    let n = values.len(); 

    let logn = (n as f32).log2() as u32;

    //let domain = EvaluationDomain::<Fb>::new(2, logn);
    let quotient_poly_degree = (2 - 1) as u64;

    // n = 2^k
    let domain_size = 1u64 << logn;    
    let mut extended_k = logn;

    while (1 << extended_k) < (domain_size * quotient_poly_degree) {
        extended_k += 1;
    }

    let mut extended_omega = Fb::root_of_unity();

    // Get extended_omega, the 2^{extended_k}'th root of unity
    // The loop computes extended_omega = omega^{2 ^ (S - extended_k)}
    // Notice that extended_omega ^ {2 ^ extended_k} = omega ^ {2^S} = 1.
    for _ in extended_k..Fb::S {
        extended_omega = extended_omega.square();
    }
    let extended_omega = extended_omega;

    // Get omega, the 2^{k}'th root of unity (i.e. n'th root of unity)
    // The loop computes omega = extended_omega ^ {2 ^ (extended_k - k)}
    //           = (omega^{2 ^ (S - extended_k)})  ^ {2 ^ (extended_k - k)}
    //           = omega ^ {2 ^ (S - k)}.
    // Notice that omega ^ {2^k} = omega ^ {2^S} = 1.

    let mut omega = extended_omega;
    for _ in logn..extended_k {
        omega = omega.square();
    }


    //let mut invs = Vec::with_capacity(4096);

    let mut evaldomain = Vec::with_capacity(4096);

    let mut acc = Fb::zero();

    let mut omega_i = Fb::one();

    let mut x_n = Fb::one();
    
    for i in 0..n {
        evaldomain.push(omega_i);

        let inv_i = (x - omega_i).invert().unwrap();

        //invs.push(inv);
        let acc_i = (values[i]) * omega_i * inv_i;

        acc += acc_i;

        omega_i = omega_i * omega;

        x_n = x_n * x;

        //println!("outside acci value{:?}",ctx.scalar_integer_ctx.assign_w(&field_to_bn(&acc_i)));
        print!("outside acc value{}:{:?} \n", i, ctx.scalar_integer_ctx.assign_w(&field_to_bn(&acc)));
    }

    acc = (x_n - Fb::one()) * Fb::from(n as u64).invert().unwrap() * acc;

    let res_expect = ctx.scalar_integer_ctx.assign_w(&field_to_bn(&acc));

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

    print!("debug: res{:?} \n", res);

    print!("debug: expect{:?} \n", res_expect);

    ctx.scalar_integer_ctx.assert_int_equal(&res, &res_expect);

    end_timer!(timer);

    run_circuit_on_bn256(ctx.into(), 22);
}