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

""" Unit tests for the unit/dimension system."""

import unittest

from steps import interface

from steps.utils import Units, Parameter

class UnitsUsage(unittest.TestCase):
    """Test units declaration and conversion."""
    def setUp(self):
        pass

    def testDeclareUnits(self):

        # Invalid input
        invInput = [
            42,
            None,
        ]
        for inp in invInput:
            with self.assertRaises(TypeError):
                Units(inp)

        # Generic tests
        for ua in Units._SI_UNITS:
            self.assertEqual(Units(ua), Units(f'({ua}^-1)^-1'))
            self.assertEqual(Units(f'm{ua}'), Units(f'(m{ua}^-1)^-1'))
            self.assertEqual(Units(f'{ua}^0'), Units(''))
            self.assertEqual(Units(f'm{ua}^2'), Units(f'u({ua}^2)'))
            self.assertEqual(Units(f'k{ua}^2'), Units(f'M({ua}^2)'))

        # Valid equivalences
        valid = [
            ('', 'm m^-1'), # Dimensionless unit
            ('m s^-1', 'mm ms^-1'),
            ('mol', 'M dm^3'),
            ('mol L^-1', 'M'),
            ('ohm', 'S^-1'),
            ('V', 'ohm A'),
            ('A', 'V ohm^-1'),
            ('A', 'V S'),
            ('F', 'A s V^-1'),
            ('kL', 'm^3'),
            ('k(g)', 'kg'),
            ('M(mg ks^-1)', 'g s^-1'),
            ('(mol m^-2)^-1 s^-1', 'mol^-1 m^2 s^-1'),
            ('m  s', 'm s'),
            ('((kg)(ms))', 'kg ms'),
            ('g^-1s', 's g^-1'),
        ]
        for val, equiv in valid:
            self.assertEqual(Units(val), Units(equiv))

        # Invalid strings
        invalid = [
            'u',
            'm/s',
            'kgs',
            'g*-1',
            'm*s^-1',
            'g**-1',
            'mols',
            '^2',
            'unknown',
            '(ms^-1',
            'm^a',
            'm^-',
            'm^',
            'Ls^-1',
            'kkg',
            'Kg',
            'm^1/2',
            'm^0.5',
        ]
        for inv in invalid:
            with self.assertRaises(Exception):
                Units(inv)

        # Invalid equivalences
        invEquivs = [
            ('m s^-1', 'm s'),
            ('M', 'mol m^-3'),
            ('m', 'mm'),
        ]
        for val, equiv in invEquivs:
            self.assertNotEqual(Units(val), Units(equiv))

        u1 = Units('s^-1')
        u1bis = Units('s^-1')
        u2 = Units('m s^-1 m^-1')
        u3 = Units('m s^-1 mm^-1')

        self.assertEqual(u1, u2)
        self.assertNotEqual(u1, u3)
        self.assertNotEqual(u2, u3)
        self.assertNotEqual(hash(u1), hash(u2))
        self.assertEqual(hash(u1), hash(u1bis))

        u1 = Units('')
        u2 = Units('m m^-1')
        u3 = Units('m mm^-1')
        self.assertTrue(u1._isDimensionless())
        self.assertTrue(u2._isDimensionless())
        self.assertFalse(u3._isDimensionless())

    def testConvertUnits(self):

        # Valid conversions
        data = [
            ('{unit}', 'm{unit}', 1e3),
            ('m{unit}', '{unit}', 1e-3),
            ('u{unit}', 'm{unit}', 1e-3),
            ('u{unit}^2 s^-1', '{unit}^2 s^-1', 1e-12),
            ('u{unit} s^-1', '{unit} s^-1', 1e-6),
            ('u{unit}^-1 s^-1', '{unit}^-1 s^-1', 1e6),
            ('u{unit}^-2 s^-1', '{unit}^-2 s^-1', 1e12),
            ('p(u{unit}^-2 s^-1)', '{unit}^-2 s^-1', 1e0),
            ('(mol m^-2)^-1 s^-1', 'mol^-1 m^2 s^-1', 1e0)
        ]

        for ua in Units._SI_UNITS:
            for src, dst, scale in data:
                u1 = Units(src.format(unit=ua))
                u2 = Units(dst.format(unit=ua))
                p1 = Parameter(1, u1, name='')
                self.assertEqual(p1.valueIn(u2), scale)

        # Invalid conversions
        data = [
            ('m', 's'),
            ('m', ''),
            ('', 'm'),
            ('L', 'kg'),
        ]
        for src, dst in data:
            u1 = Units(src)
            u2 = Units(dst)
            with self.assertRaises(Exception):
                p1 = Parameter(1, u1, name='')
                p1.valueIn(u2)

        # Incorrect types
        unit = Units('m s^-1')
        data = [
            42,
            None,
        ]
        for val in data:
            with self.assertRaises(TypeError):
                p1 = Parameter(1, u1, name='')
                p1.valueIn(val)
                
    def testRepresentations(self):

        inputStr = [
            'mol  L^-1',
            '(mol m^-2 )^-1  s^-1',
            'uM^-1  s^-1  '
        ]

        allOutput = [
            (
                [
                    'mol L^-1',
                    '(mol m^-2)^-1 s^-1',
                    'uM^-1 s^-1',
                ],
                str,
            ),
            (
                [
                    'mol L⁻¹',
                    '(mol m⁻²)⁻¹ s⁻¹',
                    'μM⁻¹ s⁻¹',
                ],
                Units._toUnicode,
            ),
            (
                [
                    '\si{mol.L^{-1}}',
                    '\si{(mol.m^{-2})^{-1}.s^{-1}}',
                    '\si{\micro M^{-1}.s^{-1}}',
                ],
                Units._toLatex,
            ),
        ]

        for outputStr, func in allOutput:
            for src, dst in zip(inputStr, outputStr):
                u = Units(src)
                self.assertEqual(func(u), dst)


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(UnitsUsage, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

