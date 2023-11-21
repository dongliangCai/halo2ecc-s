use halo2_proofs::arithmetic::CurveAffine;
use halo2_proofs::arithmetic::Field;
use halo2_proofs::arithmetic::FieldExt;

use super::base_chip::BaseChipOps;
use super::ecc_chip::EccBaseIntegerChipWrapper;
use super::ecc_chip::EccScalarIntegerChipWrapper;
use super::ecc_chip::EccChipScalarOps;
use super::integer_chip::IntegerChipOps;
use crate::assign::AssignedCondition;
use crate::assign::AssignedInteger;
use crate::circuit::ecc_chip::EccChipBaseOps;
use crate::context::GeneralScalarEccContext;
use crate::pair;
use crate::utils::field_to_bn;

impl<C: CurveAffine, N: FieldExt> EccBaseIntegerChipWrapper<C::Base, N>
    for GeneralScalarEccContext<C, N>
{
    fn base_integer_chip(&mut self) -> &mut dyn IntegerChipOps<C::Base, N> {
        &mut self.base_integer_ctx
    }
}

impl<C: CurveAffine, N: FieldExt> EccScalarIntegerChipWrapper<C::ScalarExt, N>
    for GeneralScalarEccContext<C, N>
{
    fn scalar_integer_chip(&mut self) -> &mut dyn IntegerChipOps<C::ScalarExt, N> {
        &mut self.scalar_integer_ctx
    }
}

impl<C: CurveAffine, N: FieldExt> EccChipBaseOps<C, N> for GeneralScalarEccContext<C, N> {}

impl<C: CurveAffine, N: FieldExt> EccChipScalarOps<C, N> for GeneralScalarEccContext<C, N> {
    type AssignedScalar = AssignedInteger<C::Scalar, N>;

    fn decompose_scalar<const WINDOW_SIZE: usize>(
        &mut self,
        s: &Self::AssignedScalar,
    ) -> Vec<[AssignedCondition<N>; WINDOW_SIZE]> {
        let zero = N::zero();
        let one = N::one();
        let two = one + one;
        let two_inv = two.invert().unwrap();

        let s = self.scalar_integer_ctx.reduce(&s);
        let mut bits = vec![];

        for l in s.limbs_le {
            let v = field_to_bn(&l.val);
            let mut rest = l.clone();
            for j in 0..self.scalar_integer_ctx.info.limb_bits {
                let b = self.native_ctx.borrow_mut().assign_bit(v.bit(j).into());
                let v = (rest.val - b.0.val) * two_inv;
                rest = self
                    .native_ctx
                    .borrow_mut()
                    .one_line_with_last(
                        vec![pair!(&rest, -one), pair!(&b.0, one)],
                        pair!(v, two),
                        None,
                        (vec![], None),
                    )
                    .1;
                bits.push(b);
            }

            self.native_ctx.borrow_mut().assert_constant(&rest, zero);
        }

        let padding = bits.len() % WINDOW_SIZE;
        if padding != 0 {
            let zero = self.native_ctx.borrow_mut().assign_constant(zero);
            for _ in 0..padding {
                bits.push(AssignedCondition(zero));
            }
        }

        let mut res = bits.chunks(WINDOW_SIZE)
            .map(|x| Vec::from(x).try_into().unwrap())
            .collect::<Vec<_>>();

        res.reverse();

        res
    }

    fn eval(
        &mut self, 
        values: &Vec<Self::AssignedScalar>, 
        point: &Self::AssignedScalar, 
        domain: &Vec<Self::AssignedScalar>
    ) -> Self::AssignedScalar{

        let len = values.len();

        let one = self.scalar_integer_chip().assign_int_constant(C::Scalar::one());
        
        let mut acc = self.scalar_integer_chip().assign_int_constant(C::Scalar::zero());


        let mut inv_i;
        let mut vw_i;
        let mut acc_i;
        let mut x_sub_w_i;

        let mut x = self.scalar_integer_chip().int_mul(&one, point);

        for i in 0..len {

            x_sub_w_i = self.scalar_integer_chip().int_sub(&point, &domain[i]);

            inv_i = self.scalar_integer_chip().int_unsafe_invert(&x_sub_w_i);

            vw_i = self.scalar_integer_chip().int_mul(&values[i], &domain[i]);

            acc_i = self.scalar_integer_chip().int_mul(&inv_i, &vw_i);
            
            acc = self.scalar_integer_chip().int_add(&acc, &acc_i);

            x = self.scalar_integer_chip().int_mul(&x, point);            
        }

        let x_n_sub_one = self.scalar_integer_chip().int_sub(&x, &one);
        
        let n = self.scalar_integer_chip().assign_int_constant(C::Scalar::from(len as u64));

        let (c,extract_entry) = self.scalar_integer_chip().int_div(&x_n_sub_one, &n);
        
        self.scalar_integer_ctx.ctx.borrow_mut().assert_false(&c);

        let result = self.scalar_integer_chip().int_mul(&extract_entry, &acc);

        result
    }

}
