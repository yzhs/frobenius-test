#!/bin/sh

# small_composites.sh -- generate report on false positives in LaTeX
#
# Copyright 2014 by Colin Benner <colin-software@yzhs.de>
#
# This file is part of frobenius-test.
#
# frobenius-test is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# frobenius-test is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with frobenius-test.  If not, see <http://www.gnu.org/licenses/>.

cp small_composites.template plots/small_composites.tex
sed  -i -e "s/NUMBERS/$(grep -c A all_params.log)/g" \
        -e "s/UPPER_BOUND/10000/g" \
        -e "s/FALSE_POSITIVES/$(grep -c Found all_params.log)/g" \
        -e "s/B_NOT_ZERO/$(grep -c 'b = [1-9][0-9]*' all_params.log)/g" \
        plots/small_composites.tex

