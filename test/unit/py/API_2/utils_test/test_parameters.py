####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   
###

""" Unit tests for the parameters system."""

import unittest

from steps import interface

from steps.utils import ParameterizedObject, Parameter, Units


class ParameterUsage(unittest.TestCase):
    """Test parameter declaration."""
    def setUp(self):
        pass

    def testParameterAutoNaming(self):
        param1 = Parameter(1)
        self.assertEqual(param1.name, 'param1')
        self.assertEqual(param1.value, 1)
        self.assertEqual(param1.units, None)

        param2 = \
            Parameter(
                2
            )
        self.assertEqual(param2.name, 'param2')
        self.assertEqual(param2.value, 2)
        self.assertEqual(param2.units, None)

        param3 = Parameter(
            3
        )
        self.assertEqual(param3.name, 'param3')
        self.assertEqual(param3.value, 3)
        self.assertEqual(param3.units, None)

        # Do not recognize name if it is not a direct assignement
        func = lambda x:x
        param4 = func(Parameter(4))
        self.assertEqual(param4.name, None)
        self.assertEqual(param4.value, 4)
        self.assertEqual(param4.units, None)

        param5, param6 = Parameter(5), Parameter(6)
        self.assertEqual(param5.name, None)
        self.assertEqual(param5.value, 5)
        self.assertEqual(param5.units, None)
        self.assertEqual(param6.name, None)
        self.assertEqual(param6.value, 6)
        self.assertEqual(param6.units, None)

        param9 = Parameter(9)
        param7, \
            param8 \
            = \
            Parameter(
                7
            ),\
            Parameter(
            8
            )
        self.assertEqual(param7.name, None)
        self.assertEqual(param7.value, 7)
        self.assertEqual(param7.units, None)
        self.assertEqual(param8.name, None)
        self.assertEqual(param8.value, 8)
        self.assertEqual(param8.units, None)
        self.assertEqual(param9.name, 'param9')
        self.assertEqual(param9.value, 9)
        self.assertEqual(param9.units, None)

    def testParameterDeclaration(self):
        param1 = Parameter(1, 'uM')
        self.assertEqual(param1.value, 1)
        self.assertEqual(param1.units, 'uM')

        with self.assertRaises(TypeError):
            Parameter(2, 42)

        with self.assertRaises(Exception):
            Parameter(3, 'test')

        param2 = Parameter(2, name='otherName')
        self.assertEqual(param2.name, 'otherName')
        self.assertEqual(param2.value, 2)
        self.assertEqual(param2.units, None)

        param3 = Parameter(3, DOI='DOIcontent', comment='comment content')
        self.assertEqual(param3.name, 'param3')
        self.assertEqual(param3.value, 3)
        self.assertEqual(param3.units, None)
        self.assertEqual(param3.DOI, 'DOIcontent')
        self.assertEqual(param3.comment, 'comment content')

        with self.assertRaises(AttributeError):
            param3.unknownAttr

        with self.assertRaises(AttributeError):
            param3.unknownAttr = 42

        with self.assertRaises(Exception):
            param3.comment = 'other comment'

        with self.assertRaises(Exception):
            Parameter(4, value='issue')

        with self.assertRaises(TypeError):
            Parameter(5, name=42)

        with self.assertRaises(Exception):
            Parameter(5, _value=10)

        param1 = Parameter(5, 'M')
        param2 = Parameter(10, 's')
        param3 = Parameter(param1 / param2)
        param4 = Parameter(param1 / param2, 'uM s^-1')

        self.assertEqual(param3.name, 'param3')
        self.assertEqual(param3.value, 0.5)
        self.assertEqual(param3.units, 'M s^-1')
        self.assertEqual(param3._fullname, 'param1 / param2')

        self.assertEqual(param4.name, 'param4')
        self.assertEqual(param4.value, 500000)
        self.assertEqual(param4.units, 'uM s^-1')
        self.assertEqual(param4._fullname, 'param1 / param2')


    def testParameterConversions(self):
        param1 = Parameter(1, 'uM')
        self.assertEqual(param1.valueIn(Units('M')), 1e-6)
        self.assertEqual(param1._valueInSI(), 1e-3)

        with self.assertRaises(Exception):
            param1.valueIn(Units('m'))

        with self.assertRaises(Exception):
            param1.valueIn(Units('M^-1'))

        with self.assertRaises(TypeError):
            param1.valueIn(42)

        param2 = Parameter(2, 'uM^-1')
        self.assertEqual(param2.valueIn(Units('M^-1')), 2e6)
        self.assertEqual(param2._valueInSI(), 2e3)

        param3 = Parameter(3)
        with self.assertRaises(Exception):
            param2.valueIn(Units('m'))
        self.assertEqual(param3.valueIn(None), 3)
        self.assertEqual(param3._valueInSI(), 3)

        param1 = Parameter('test')
        with self.assertRaises(Exception):
            param1.valueIn('s')
        with self.assertRaises(Exception):
            param1.valueInSI()

    def testParameterArithmetic(self):

        # CHeck that autonaming is not accidentally picking up 'val' as name for Parameter(10)
        val = Parameter(10, '') + Parameter(11, '')
        self.assertFalse(val._isNamed())


        # With numbers
        param2 = Parameter(2, 'uM')
        param5 = 5

        add = param2 + param5
        sub = param2 - param5
        mul = param2 * param5
        div = param2 / param5

        add_inv = param5 + param2
        sub_inv = param5 - param2
        mul_inv = param5 * param2
        div_inv = param5 / param2

        power = param2 ** 2

        self.assertTrue(param2._isUserDefined())
        self.assertFalse(add._isUserDefined())

        self.assertEqual(add.value, 7)
        self.assertEqual(add.units, 'uM')
        self.assertEqual(add.name, 'param2 + 5')

        self.assertEqual(sub.value, -3)
        self.assertEqual(sub.units, 'uM')
        self.assertEqual(sub.name, 'param2 - 5')

        self.assertEqual(mul.value, 10)
        self.assertEqual(mul.units, 'uM')
        self.assertEqual(mul.name, 'param2 * 5')

        self.assertEqual(div.value, 0.4)
        self.assertEqual(div.units, 'uM')
        self.assertEqual(div.name, 'param2 / 5')

        self.assertEqual(add_inv.value, 7)
        self.assertEqual(add_inv.units, 'uM')
        self.assertEqual(add_inv.name, '5 + param2')

        self.assertEqual(sub_inv.value, 3)
        self.assertEqual(sub_inv.units, 'uM')
        self.assertEqual(sub_inv.name, '5 - param2')

        self.assertEqual(mul_inv.value, 10)
        self.assertEqual(mul_inv.units, 'uM')
        self.assertEqual(mul_inv.name, '5 * param2')

        self.assertEqual(div_inv.value, 2.5)
        self.assertEqual(div_inv.units, 'uM^-1')
        self.assertEqual(div_inv.name, '5 / param2')

        self.assertEqual(power.value, 4)
        self.assertEqual(power.units, 'uM^2')
        self.assertEqual(power.name, 'param2 ** 2')

        # No units
        param2 = Parameter(2, '')
        param5 = Parameter(5, '')

        add = param2 + param5
        sub = param2 - param5
        mul = param2 * param5
        div = param2 / param5
        power = param2 ** 2

        self.assertEqual(add.value, 7)
        self.assertEqual(add.units, '')
        self.assertEqual(add.name, 'param2 + param5')

        self.assertEqual(sub.value, -3)
        self.assertEqual(sub.units, '')
        self.assertEqual(sub.name, 'param2 - param5')

        self.assertEqual(mul.value, 10)
        self.assertEqual(mul.units, '')
        self.assertEqual(mul.name, 'param2 * param5')

        self.assertEqual(div.value, 0.4)
        self.assertEqual(div.units, '')
        self.assertEqual(div.name, 'param2 / param5')

        self.assertEqual(power.value, 4)
        self.assertEqual(power.units, '')
        self.assertEqual(power.name, 'param2 ** 2')

        # Forbidden operations
        with self.assertRaises(Exception):
            param2 * Parameter(10)

        # Units already raised to some power
        param2 = Parameter(2, 'uM^2')

        power = param2 ** -2

        self.assertEqual(power.value, 0.25)
        self.assertEqual(power.units, 'uM^-4')
        self.assertEqual(power.name, 'param2 ** -2')

        # Two named Parameter objects, incompatible units
        param2 = Parameter(2, 'm s^-1')
        param5 = Parameter(5, 'uM')

        with self.assertRaises(Exception):
            add = param2 + param5
        with self.assertRaises(Exception):
            invadd = param5 + param2
        with self.assertRaises(Exception):
            sub = param2 - param5
        with self.assertRaises(Exception):
            invsub = param5 - param2

        mul = param2 * param5
        invmul = param5 * param2
        div = param2 / param5
        invdiv = param5 / param2

        power = param2 ** 2

        self.assertEqual(mul.value, 10)
        self.assertEqual(mul.units, 'm s^-1 uM')
        self.assertEqual(mul.name, 'param2 * param5')
        self.assertEqual(invmul.value, 10)
        self.assertEqual(invmul.units, 'uM m s^-1')
        self.assertEqual(invmul.name, 'param5 * param2')

        self.assertEqual(div.value, 0.4)
        self.assertEqual(div.units, 'm s^-1 uM^-1')
        self.assertEqual(div.name, 'param2 / param5')
        self.assertEqual(invdiv.value, 2.5)
        self.assertEqual(invdiv.units, 'uM (m s^-1)^-1')
        self.assertEqual(invdiv.name, 'param5 / param2')

        self.assertEqual(power.value, 4)
        self.assertEqual(power.units, '(m s^-1)^2')
        self.assertEqual(power.name, 'param2 ** 2')

        # Two named Parameter objects, compatible units
        param2 = Parameter(2, 'm s^-1')
        param5 = Parameter(5, 'um ms^-1')

        add = param2 + param5
        invadd = param5 + param2
        sub = param2 - param5
        invsub = param5 - param2

        mul = param2 * param5
        invmul = param5 * param2
        div = param2 / param5
        invdiv = param5 / param2

        self.assertEqual(add.value, 2 + 5e-3)
        self.assertEqual(add.units, 'm s^-1')
        self.assertEqual(add.name, 'param2 + param5')
        self.assertEqual(invadd.value, 2005)
        self.assertEqual(invadd.units, 'um ms^-1')
        self.assertEqual(invadd.name, 'param5 + param2')

        self.assertEqual(sub.value, 2 - 5e-3)
        self.assertEqual(sub.units, 'm s^-1')
        self.assertEqual(sub.name, 'param2 - param5')
        self.assertEqual(invsub.value, -1995)
        self.assertEqual(invsub.units, 'um ms^-1')
        self.assertEqual(invsub.name, 'param5 - param2')

        self.assertEqual(mul.value, 10)
        self.assertEqual(mul.units, 'm s^-1 um ms^-1')
        self.assertEqual(mul.name, 'param2 * param5')
        self.assertEqual(invmul.value, 10)
        self.assertEqual(invmul.units, 'um ms^-1 m s^-1')
        self.assertEqual(invmul.name, 'param5 * param2')

        self.assertEqual(div.value, 0.4)
        self.assertEqual(div.units, 'm s^-1 (um ms^-1)^-1')
        self.assertEqual(div.name, 'param2 / param5')
        self.assertEqual(invdiv.value, 2.5)
        self.assertEqual(invdiv.units, 'um ms^-1 (m s^-1)^-1')
        self.assertEqual(invdiv.name, 'param5 / param2')

        # Two Parameter objects, one named, incompatible units
        param2 = Parameter(2, 'm s^-1')

        with self.assertRaises(Exception):
            add = param2 + Parameter(5, 'uM')
        with self.assertRaises(Exception):
            invadd = Parameter(5, 'uM') + param2
        with self.assertRaises(Exception):
            sub = param2 - Parameter(5, 'uM')
        with self.assertRaises(Exception):
            invsub = Parameter(5, 'uM') - param2

        mul = param2 * Parameter(5, 'uM')
        invmul = Parameter(5, 'uM') * param2
        div = param2 / Parameter(5, 'uM')
        invdiv = Parameter(5, 'uM') / param2

        self.assertEqual(mul.value, 10)
        self.assertEqual(mul.units, 'm s^-1 uM')
        self.assertEqual(mul.name, 'param2 * (5 uM)')
        self.assertEqual(invmul.value, 10)
        self.assertEqual(invmul.units, 'uM m s^-1')
        self.assertEqual(invmul.name, '(5 uM) * param2')

        self.assertEqual(div.value, 0.4)
        self.assertEqual(div.units, 'm s^-1 uM^-1')
        self.assertEqual(div.name, 'param2 / (5 uM)')
        self.assertEqual(invdiv.value, 2.5)
        self.assertEqual(invdiv.units, 'uM (m s^-1)^-1')
        self.assertEqual(invdiv.name, '(5 uM) / param2')

        # Two named Parameter objects, compatible units
        param2 = Parameter(2, 'm s^-1')

        add = param2 + Parameter(5, 'um ms^-1')
        invadd = Parameter(5, 'um ms^-1') + param2
        sub = param2 - Parameter(5, 'um ms^-1')
        invsub = Parameter(5, 'um ms^-1') - param2

        mul = param2 * Parameter(5, 'um ms^-1')
        invmul = Parameter(5, 'um ms^-1') * param2
        div = param2 / Parameter(5, 'um ms^-1')
        invdiv = Parameter(5, 'um ms^-1') / param2

        self.assertEqual(add.value, 2 + 5e-3)
        self.assertEqual(add.units, 'm s^-1')
        self.assertEqual(add.name, 'param2 + 0.005')
        self.assertEqual(invadd.value, 5e-3 + 2)
        self.assertEqual(invadd.units, 'm s^-1')
        self.assertEqual(invadd.name, '0.005 + param2')

        self.assertEqual(sub.value, 2 - 5e-3)
        self.assertEqual(sub.units, 'm s^-1')
        self.assertEqual(sub.name, 'param2 - 0.005')
        self.assertEqual(invsub.value, 5e-3 - 2)
        self.assertEqual(invsub.units, 'm s^-1')
        self.assertEqual(invsub.name, '0.005 - param2')

        self.assertEqual(mul.value, 10)
        self.assertEqual(mul.units, 'm s^-1 um ms^-1')
        self.assertEqual(mul.name, 'param2 * (5 um ms^-1)')
        self.assertEqual(invmul.value, 10)
        self.assertEqual(invmul.units, 'um ms^-1 m s^-1')
        self.assertEqual(invmul.name, '(5 um ms^-1) * param2')

        self.assertEqual(div.value, 0.4)
        self.assertEqual(div.units, 'm s^-1 (um ms^-1)^-1')
        self.assertEqual(div.name, 'param2 / (5 um ms^-1)')
        self.assertEqual(invdiv.value, 2.5)
        self.assertEqual(invdiv.units, 'um ms^-1 (m s^-1)^-1')
        self.assertEqual(invdiv.name, '(5 um ms^-1) / param2')

        # Two unnamed Parameter objects, incompatible units
        with self.assertRaises(Exception):
            add = Parameter(2, 'm s^-1') + Parameter(5, 'uM')
        with self.assertRaises(Exception):
            invadd = Parameter(5, 'uM') + Parameter(2, 'm s^-1')
        with self.assertRaises(Exception):
            sub = Parameter(2, 'm s^-1') - Parameter(5, 'uM')
        with self.assertRaises(Exception):
            invsub = Parameter(5, 'uM') - Parameter(2, 'm s^-1')

        mul = Parameter(2, 'm s^-1') * Parameter(5, 'uM')
        invmul = Parameter(5, 'uM') * Parameter(2, 'm s^-1')
        div = Parameter(2, 'm s^-1') / Parameter(5, 'uM')
        invdiv = Parameter(5, 'uM') / Parameter(2, 'm s^-1')

        power = Parameter(2, 'm s^-1') ** 2

        self.assertEqual(mul.value, 10)
        self.assertEqual(mul.units, 'm s^-1 uM')
        self.assertEqual(mul.name, None)
        self.assertEqual(invmul.value, 10)
        self.assertEqual(invmul.units, 'uM m s^-1')
        self.assertEqual(invmul.name, None)

        self.assertEqual(div.value, 0.4)
        self.assertEqual(div.units, 'm s^-1 uM^-1')
        self.assertEqual(div.name, None)
        self.assertEqual(invdiv.value, 2.5)
        self.assertEqual(invdiv.units, 'uM (m s^-1)^-1')
        self.assertEqual(invdiv.name, None)

        self.assertEqual(power.value, 4)
        self.assertEqual(power.units, '(m s^-1)^2')
        self.assertEqual(power.name, None)

        # Two unnamed Parameter objects, compatible units
        add = Parameter(2, 'm s^-1') + Parameter(5, 'um ms^-1')
        invadd = Parameter(5, 'um ms^-1') + Parameter(2, 'm s^-1')
        sub = Parameter(2, 'm s^-1') - Parameter(5, 'um ms^-1')
        invsub = Parameter(5, 'um ms^-1') - Parameter(2, 'm s^-1')

        mul = Parameter(2, 'm s^-1') * Parameter(5, 'um ms^-1')
        invmul = Parameter(5, 'um ms^-1') * Parameter(2, 'm s^-1')
        div = Parameter(2, 'm s^-1') / Parameter(5, 'um ms^-1')
        invdiv = Parameter(5, 'um ms^-1') / Parameter(2, 'm s^-1')

        self.assertEqual(add.value, 2 + 5e-3)
        self.assertEqual(add.units, 'm s^-1')
        self.assertEqual(add.name, None)
        self.assertEqual(invadd.value, 5 + 2e3)
        self.assertEqual(invadd.units, 'um ms^-1')
        self.assertEqual(invadd.name, None)

        self.assertEqual(sub.value, 2 - 5e-3)
        self.assertEqual(sub.units, 'm s^-1')
        self.assertEqual(sub.name, None)
        self.assertEqual(invsub.value, 5 - 2e3)
        self.assertEqual(invsub.units, 'um ms^-1')
        self.assertEqual(invsub.name, None)

        self.assertEqual(mul.value, 10)
        self.assertEqual(mul.units, 'm s^-1 um ms^-1')
        self.assertEqual(mul.name, None)
        self.assertEqual(invmul.value, 10)
        self.assertEqual(invmul.units, 'um ms^-1 m s^-1')
        self.assertEqual(invmul.name, None)

        self.assertEqual(div.value, 0.4)
        self.assertEqual(div.units, 'm s^-1 (um ms^-1)^-1')
        self.assertEqual(div.name, None)
        self.assertEqual(invdiv.value, 2.5)
        self.assertEqual(invdiv.units, 'um ms^-1 (m s^-1)^-1')
        self.assertEqual(invdiv.name, None)

        # More than one operation, all named, compatible units
        param2 = Parameter(2, 'm s^-1')
        param4 = Parameter(4, 'cm ds^-1')
        param5 = Parameter(5, 'um ms^-1')

        add_add = param2 + param4 + param5
        add_add_inv = param5 + param4 + param2

        sub_sub = param2 - param4 - param5
        sub_sub_inv = param5 - param4 - param2

        mul_mul = param2 * param4 * param5
        mul_mul_inv = param5 * param4 * param2

        div_div = param2 / param4 / param5
        div_div_inv = param5 / param4 / param2

        add_sub = param2 + param4 - param5
        add_sub_inv = param5 + param4 - param2

        mul_div = param2 * param4 / param5
        mul_div_inv = param5 * param4 / param2

        add_power = (param2 + param5) ** 2
        sub_power = (param2 - param5) ** 2
        mul_power = (param2 * param5) ** 2
        div_power = (param2 / param5) ** 2

        self.assertEqual(add_add.value, 2 + 4e-1 + 5e-3)
        self.assertEqual(add_add.units, 'm s^-1')
        self.assertEqual(add_add.name, 'param2 + param4 + param5')
        self.assertEqual(add_add_inv.value, 2e3 + 4e2 + 5)
        self.assertEqual(add_add_inv.units, 'um ms^-1')
        self.assertEqual(add_add_inv.name, 'param5 + param4 + param2')

        self.assertEqual(sub_sub.value, 2 - 4e-1 - 5e-3)
        self.assertEqual(sub_sub.units, 'm s^-1')
        self.assertEqual(sub_sub.name, 'param2 - param4 - param5')
        self.assertEqual(sub_sub_inv.value, 5 - 4e2 - 2e3)
        self.assertEqual(sub_sub_inv.units, 'um ms^-1')
        self.assertEqual(sub_sub_inv.name, 'param5 - param4 - param2')

        self.assertEqual(mul_mul.value, 40)
        self.assertEqual(mul_mul.units, 'm s^-1 cm ds^-1 um ms^-1')
        self.assertEqual(mul_mul.name, 'param2 * param4 * param5')
        self.assertEqual(mul_mul_inv.value, 40)
        self.assertEqual(mul_mul_inv.units, 'um ms^-1 cm ds^-1 m s^-1')
        self.assertEqual(mul_mul_inv.name, 'param5 * param4 * param2')

        self.assertEqual(div_div.value, 0.1)
        self.assertEqual(div_div.units, 'm s^-1 (cm ds^-1)^-1 (um ms^-1)^-1')
        self.assertEqual(div_div.name, 'param2 / param4 / param5')
        self.assertEqual(div_div_inv.value, 0.625)
        self.assertEqual(div_div_inv.units, 'um ms^-1 (cm ds^-1)^-1 (m s^-1)^-1')
        self.assertEqual(div_div_inv.name, 'param5 / param4 / param2')

        self.assertEqual(add_power.value, (2 + 5e-3) ** 2)
        self.assertEqual(add_power.units, '(m s^-1)^2')
        self.assertEqual(add_power.name, '(param2 + param5) ** 2')

        self.assertEqual(sub_power.value, (2 - 5e-3) ** 2)
        self.assertEqual(sub_power.units, '(m s^-1)^2')
        self.assertEqual(sub_power.name, '(param2 - param5) ** 2')

        self.assertEqual(mul_power.value, 10 ** 2)
        self.assertEqual(mul_power.units, '(m s^-1 um ms^-1)^2')
        self.assertEqual(mul_power.name, '(param2 * param5) ** 2')

        self.assertEqual(div_power.value, 0.4 ** 2)
        self.assertEqual(div_power.units, '(m s^-1 (um ms^-1)^-1)^2')
        self.assertEqual(div_power.name, '(param2 / param5) ** 2')

        param2 = Parameter(2, 'm s^-1')
        param5a = Parameter(5, 'um')
        param5b = Parameter(1, 'ms^-1')
        param5c = Parameter(1, 'ms')

        add_mul = param2 + param5a * param5b
        add_mul_inv = param5a * param5b + param2

        add_div = param2 + param5a / param5c
        add_div_inv = param5a / param5c + param2

        sub_mul = param2 - param5a * param5b
        sub_mul_inv = param5a * param5b - param2

        sub_div = param2 - param5a / param5c
        sub_div_inv = param5a / param5c - param2

        self.assertEqual(add_mul.value, 2 + 5e-3)
        self.assertEqual(add_mul.units, 'm s^-1')
        self.assertEqual(add_mul.name, 'param2 + param5a * param5b')
        self.assertEqual(add_mul_inv.value, 2e3 + 5)
        self.assertEqual(add_mul_inv.units, 'um ms^-1')
        self.assertEqual(add_mul_inv.name, 'param5a * param5b + param2')

        self.assertEqual(add_div.value, 2 + 5e-3)
        self.assertEqual(add_div.units, 'm s^-1')
        self.assertEqual(add_div.name, 'param2 + param5a / param5c')
        self.assertEqual(add_div_inv.value, 2e3 + 5)
        self.assertEqual(add_div_inv.units, 'um ms^-1')
        self.assertEqual(add_div_inv.name, 'param5a / param5c + param2')

        self.assertEqual(sub_mul.value, 2 - 5e-3)
        self.assertEqual(sub_mul.units, 'm s^-1')
        self.assertEqual(sub_mul.name, 'param2 - param5a * param5b')
        self.assertEqual(sub_mul_inv.value, -1995)
        self.assertEqual(sub_mul_inv.units, 'um ms^-1')
        self.assertEqual(sub_mul_inv.name, 'param5a * param5b - param2')

        self.assertEqual(sub_div.value, 2 - 5e-3)
        self.assertEqual(sub_div.units, 'm s^-1')
        self.assertEqual(sub_div.name, 'param2 - param5a / param5c')
        self.assertEqual(sub_div_inv.value, -1995)
        self.assertEqual(sub_div_inv.units, 'um ms^-1')
        self.assertEqual(sub_div_inv.name, 'param5a / param5c - param2')

        #Parentheses
        param2a = Parameter(1, 'm s^-1')
        param2b = Parameter(1, 'mm ms^-1')
        param2c = Parameter(3, 'm s^-1')
        param5 = Parameter(5, 'um ms^-1')

        mul_add = (param2a + param2b) * param5
        mul_add_inv = param5 * (param2a + param2b)

        div_add = (param2a + param2b) / param5
        div_add_inv = param5 / (param2a + param2b)

        mul_sub = (param2c - param2b) * param5
        mul_sub_inv = param5 * (param2c - param2b)

        div_sub = (param2c - param2b) / param5
        div_sub_inv = param5 / (param2c - param2b)

        self.assertEqual(mul_add.value, 10)
        self.assertEqual(mul_add.units, 'm s^-1 um ms^-1')
        self.assertEqual(mul_add.name, '(param2a + param2b) * param5')
        self.assertEqual(mul_add_inv.value, 10)
        self.assertEqual(mul_add_inv.units, 'um ms^-1 m s^-1')
        self.assertEqual(mul_add_inv.name, 'param5 * (param2a + param2b)')

        self.assertEqual(div_add.value, 0.4)
        self.assertEqual(div_add.units, 'm s^-1 (um ms^-1)^-1')
        self.assertEqual(div_add.name, '(param2a + param2b) / param5')
        self.assertEqual(div_add_inv.value, 2.5)
        self.assertEqual(div_add_inv.units, 'um ms^-1 (m s^-1)^-1')
        self.assertEqual(div_add_inv.name, 'param5 / (param2a + param2b)')

        self.assertEqual(mul_sub.value, 10)
        self.assertEqual(mul_sub.units, 'm s^-1 um ms^-1')
        self.assertEqual(mul_sub.name, '(param2c - param2b) * param5')
        self.assertEqual(mul_sub_inv.value, 10)
        self.assertEqual(mul_sub_inv.units, 'um ms^-1 m s^-1')
        self.assertEqual(mul_sub_inv.name, 'param5 * (param2c - param2b)')

        self.assertEqual(div_sub.value, 0.4)
        self.assertEqual(div_sub.units, 'm s^-1 (um ms^-1)^-1')
        self.assertEqual(div_sub.name, '(param2c - param2b) / param5')
        self.assertEqual(div_sub_inv.value, 2.5)
        self.assertEqual(div_sub_inv.units, 'um ms^-1 (m s^-1)^-1')
        self.assertEqual(div_sub_inv.name, 'param5 / (param2c - param2b)')

        self.assertEqual(
            set(mul_add._getAllDependencies()),
            set([
                param2a + param2b,
                param2a,
                param2b,
                param5,
            ]),
        )

        self.assertEqual(
            set(mul_sub_inv._getAllDependencies()),
            set([
                param2c - param2b,
                param2b,
                param2c,
                param5,
            ]),
        )

        composed = mul_add * mul_sub_inv
        self.assertEqual(
            set(composed._getAllDependencies()),
            set([
                mul_add,
                mul_sub_inv,
                param2a + param2b,
                param2c - param2b,
                param2a,
                param2b,
                param2c,
                param5,
            ]),
        )

        with self.assertRaises(TypeError):
            param5 + 'test'

        with self.assertRaises(TypeError):
            param5 - 'test'

        with self.assertRaises(TypeError):
            param5 * 'test'

        with self.assertRaises(TypeError):
            param5 / 'test'

        with self.assertRaises(TypeError):
            param5 ** 0.5

        with self.assertRaises(TypeError):
            param5 ** 'test'


    def testParameterizedDeclaration(self):
        class TestClass(ParameterizedObject):
            def __init__(self, a, b, c, setb=True, units=Units(''), *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.a = a
                if setb:
                    self.b = b
                self._units = units
                self.c = c
                self.d = None

            def getUnits(self):
                return self._units

            @property
            @ParameterizedObject.RegisterGetter(Units('m s^-1'))
            def a(self):
                self.d = 'geta'
                return 42

            @a.setter
            @ParameterizedObject.RegisterSetter(Units('m s^-1'))
            def a(self, v):
                self.d = 'seta'
                self.v = v

            @property
            @ParameterizedObject.RegisterGetter(Units('mV m^-1'))
            def b(self):
                self.d = 'getb'
                return 21

            @b.setter
            @ParameterizedObject.RegisterSetter(Units('mV m^-1'))
            def b(self, v):
                self.d = 'setb'
                self.v = v

            @property
            @ParameterizedObject.RegisterGetter(getUnits)
            def c(self):
                self.d = 'getc'
                return 21

            @c.setter
            @ParameterizedObject.RegisterSetter(getUnits)
            def c(self, v):
                self.d = 'setc'
                self.v = v

        TestClass.RegisterParameter('e', Units('S'), 123)

        po = TestClass(1, 2, 3, True)
        self.assertEqual(po.a, 1)
        self.assertEqual(po.b, 2)
        self.assertEqual(po.d, None)
        po.a = 5
        self.assertEqual(po.d, 'seta')
        po.b = 10
        self.assertEqual(po.d, 'setb')
        po.d = None
        self.assertEqual(po.a, 5)
        self.assertEqual(po.b, 10)
        self.assertEqual(po.d, None)

        po = TestClass(1, None, None, False)
        self.assertEqual(po.a, 1)
        self.assertEqual(po.d, None)
        self.assertEqual(po.b, 21)
        self.assertEqual(po.d, 'getb')
        po.b = 10
        self.assertEqual(po.d, 'setb')
        self.assertEqual(po.b, 10)

        param1  = Parameter(1, 'mm s^-1')
        param2  = Parameter(1, 'V mm^-1')

        param1b = Parameter(1, 'm s^-1')
        param2b = Parameter(1, 'mV m^-1')

        param3 = Parameter(1, 'M')

        po = TestClass(param1, param2, param3, True, Units('mol L^-1'))
        self.assertEqual(po.a, 1e-3)
        self.assertEqual(po.b, 1e6)
        self.assertEqual(po.c, 1)
        self.assertEqual(po.d, None)
        po.a = param1b
        self.assertEqual(po.a, 1)
        self.assertEqual(po.d, 'seta')
        po.b = param2b
        self.assertEqual(po.b, 1)
        self.assertEqual(po.d, 'setb')
        po.d = None
        with self.assertRaises(Exception):
            po.a = param2
        self.assertEqual(po.d, None)
        with self.assertRaises(Exception):
            po.b = param1
        self.assertEqual(po.d, None)

        with self.assertRaises(Exception):
            TestClass(param1, param2, param3, True, Units('V'))

        param1 = Parameter(1)
        po = TestClass(param1, None, None, False)
        self.assertEqual(po.a, 1)
        self.assertEqual(param1.units, 'm s^-1')
        with self.assertRaises(Exception):
            po.b = param1
        self.assertEqual(po.d, None)
        self.assertEqual(po.b, 21)
        self.assertEqual(po.d, 'getb')

        # e was not set yet but it should return its default value
        self.assertEqual(po.e, 123)
        po.e = Parameter(1, 'mS')
        self.assertEqual(po.e, 1e-3)
        with self.assertRaises(Exception):
            po.e = Parameter(1, 'V')

        param4 = Parameter(po)
        self.assertEqual(param4.valueIn('m'), po)

        param2 = Parameter(2, 'm s^-1')
        po._setParameter('param2', param2)
        po._setParameter('param3', 3)

        param4 = Parameter(4, 'M V^-1 s^-1')
        po._setAdvancedParameter((('col10', 10), ('col20', 20)), param4)
        po._setAdvancedParameter((('col10', 100), ('col20', 200)), 5)

        self.assertEqual(
            set(po._getAllParams()),
            set([
                Parameter(None, '', name=''),
                param1,
                Parameter(1, 'mS'),
                param2,
                Parameter(3, name=''),
                param4,
                Parameter(5, name=''),
            ])
        )

    def testExplicitCasts(self):
        param1 = Parameter(1.3e3, 'uM')
        param2 = Parameter(2454.5, 'ms')
        param3 = Parameter(3.8e3, 'uM')
        param4 = Parameter(4784.5, 'ms')

        self.assertEqual(int(param1), 1)
        self.assertEqual(round(param1), 1.0)
        self.assertAlmostEqual(float(param1), 1.3)

        self.assertEqual(int(param2), 2)
        self.assertEqual(round(param2), 2.0)
        self.assertAlmostEqual(float(param2), 2.4545)

        self.assertEqual(int(param3), 3)
        self.assertEqual(round(param3), 4.0)
        self.assertAlmostEqual(float(param3), 3.8)

        self.assertEqual(int(param4), 4)
        self.assertEqual(round(param4), 5.0)
        self.assertAlmostEqual(float(param4), 4.7845)

        param6 = Parameter(True)
        param7 = Parameter(False)

        self.assertTrue(bool(param6))
        self.assertFalse(bool(param7))



def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(ParameterUsage, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

