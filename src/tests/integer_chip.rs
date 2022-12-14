use crate::circuit::base_chip::{BaseChip, BaseChipConfig, BaseChipOps};
use crate::circuit::integer_chip::IntegerChipOps;
use crate::circuit::range_chip::RangeChip;
use crate::circuit::range_chip::RangeChipConfig;
use crate::context::Records;
use crate::context::{Context, IntegerContext};
use crate::tests::{random_bls12_381_fq, random_bls12_381_fr, random_fq};
use crate::utils::field_to_bn;
use ark_std::{end_timer, start_timer};
use halo2_proofs::arithmetic::{BaseExt, Field, FieldExt};
use halo2_proofs::pairing::bn256::Fr;
use halo2_proofs::{
    circuit::{Layouter, SimpleFloorPlanner},
    dev::MockProver,
    plonk::{Circuit, ConstraintSystem, Error},
};
use std::cell::RefCell;
use std::marker::PhantomData;
use std::sync::Arc;

#[derive(Clone)]
struct TestChipConfig {
    base_chip_config: BaseChipConfig,
    range_chip_config: RangeChipConfig,
}

#[derive(Default)]
struct TestCircuit<W: BaseExt, N: FieldExt> {
    records: Records<N>,
    _phantom: PhantomData<W>,
}

impl<W: BaseExt, N: FieldExt> Circuit<N> for TestCircuit<W, N> {
    type Config = TestChipConfig;
    type FloorPlanner = SimpleFloorPlanner;

    fn without_witnesses(&self) -> Self {
        Self::default()
    }

    fn configure(meta: &mut ConstraintSystem<N>) -> Self::Config {
        let base_chip_config = BaseChip::configure(meta);
        let range_chip_config = RangeChip::<W, N>::configure(meta);
        TestChipConfig {
            base_chip_config,
            range_chip_config,
        }
    }

    fn synthesize(
        &self,
        config: Self::Config,
        mut layouter: impl Layouter<N>,
    ) -> Result<(), Error> {
        let base_chip = BaseChip::new(config.base_chip_config);
        let range_chip = RangeChip::<W, N>::new(config.range_chip_config);

        range_chip.init_table(&mut layouter)?;

        layouter.assign_region(
            || "base",
            |mut region| {
                let timer = start_timer!(|| "assign");
                self.records
                    .assign_all(&mut region, &base_chip, &range_chip)?;
                end_timer!(timer);
                Ok(())
            },
        )?;

        Ok(())
    }
}

#[test]
fn test_integer_chip_st() {
    let ctx = RefCell::new(Context::new());
    let mut ctx = IntegerContext::<halo2_proofs::pairing::bn256::Fq, Fr>::new(ctx);

    let a = random_fq();
    let mut b = random_fq();

    let mut timeout = 3;
    while b.is_zero().into() && timeout > 0 {
        b = random_fq();
        timeout -= 1;
    }
    assert!(timeout > 0);

    let c = a + b;
    let d = a - b;
    let e = a * b;
    let f = a * b.invert().unwrap();

    let a = ctx.assign_w(&field_to_bn(&a));
    let b = ctx.assign_w(&field_to_bn(&b));

    let c1 = ctx.assign_w(&field_to_bn(&c));
    let c2 = ctx.int_add(&a, &b);
    ctx.assert_int_equal(&c1, &c2);

    let d1 = ctx.assign_w(&field_to_bn(&d));
    let d2 = ctx.int_sub(&a, &b);
    ctx.assert_int_equal(&d1, &d2);

    let e1 = ctx.assign_w(&field_to_bn(&e));
    let e2 = ctx.int_mul(&a, &b);
    ctx.assert_int_equal(&e1, &e2);

    let f1 = ctx.assign_w(&field_to_bn(&f));
    let f2 = ctx.int_div(&a, &b).1;
    ctx.assert_int_equal(&f1, &f2);

    let zero = ctx.int_sub(&a, &a);
    let (g1, _) = ctx.int_div(&a, &zero);
    ctx.ctx.borrow_mut().assert_true(&g1);

    const K: u32 = 20;
    let circuit = TestCircuit::<halo2_proofs::pairing::bn256::Fq, Fr> {
        records: Arc::try_unwrap(ctx.ctx.into_inner().records)
            .unwrap()
            .into_inner()
            .unwrap(),
        _phantom: PhantomData,
    };
    let prover = match MockProver::run(K, &circuit, vec![]) {
        Ok(prover) => prover,
        Err(e) => panic!("{:#?}", e),
    };
    assert_eq!(prover.verify(), Ok(()));
}

#[test]
fn test_bls12_381_fq_integer_chip_st() {
    let ctx = RefCell::new(Context::new());
    let mut ctx = IntegerContext::<halo2_proofs::pairing::bls12_381::Fq, Fr>::new(ctx);

    for _ in 0..1000 {
        let a = random_bls12_381_fq();
        let b = random_bls12_381_fq();
        let ab = a * b;

        let a = ctx.assign_w(&field_to_bn(&a));
        let b = ctx.assign_w(&field_to_bn(&b));
        let ab0 = ctx.assign_w(&field_to_bn(&ab));

        let ab1 = ctx.int_mul(&a, &b);

        ctx.assert_int_equal(&ab0, &ab1);
    }

    const K: u32 = 20;
    let circuit = TestCircuit::<halo2_proofs::pairing::bls12_381::Fq, Fr> {
        records: Arc::try_unwrap(ctx.ctx.into_inner().records)
            .unwrap()
            .into_inner()
            .unwrap(),
        _phantom: PhantomData,
    };
    let prover = match MockProver::run(K, &circuit, vec![]) {
        Ok(prover) => prover,
        Err(e) => panic!("{:#?}", e),
    };
    assert_eq!(prover.verify(), Ok(()));
}

#[test]
fn test_bls12_381_fr_integer_chip_st() {
    let ctx = RefCell::new(Context::new());
    let mut ctx = IntegerContext::<halo2_proofs::pairing::bls12_381::Fr, Fr>::new(ctx);

    for _ in 0..1000 {
        let a = random_bls12_381_fr();
        let b = random_bls12_381_fr();
        let ab = a * b;

        let a = ctx.assign_w(&field_to_bn(&a));
        let b = ctx.assign_w(&field_to_bn(&b));
        let ab0 = ctx.assign_w(&field_to_bn(&ab));

        let ab1 = ctx.int_mul(&a, &b);

        ctx.assert_int_equal(&ab0, &ab1);
    }

    const K: u32 = 20;
    let circuit = TestCircuit::<halo2_proofs::pairing::bls12_381::Fr, Fr> {
        records: Arc::try_unwrap(ctx.ctx.into_inner().records)
            .unwrap()
            .into_inner()
            .unwrap(),
        _phantom: PhantomData,
    };
    let prover = match MockProver::run(K, &circuit, vec![]) {
        Ok(prover) => prover,
        Err(e) => panic!("{:#?}", e),
    };
    assert_eq!(prover.verify(), Ok(()));
}
