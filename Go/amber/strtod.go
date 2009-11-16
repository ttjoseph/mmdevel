package amber
// This source code was downloaded from
// http://www.jbox.dk/sanos/source/lib/strtod.c.html
// and ported to Go by Tom Joseph <thomas.joseph@mssm.edu>

// Files used:
// strtod.c: Convert string to double 
// float.h: Constants for floating point values

// Copyright (C) 2002 Michael Ringgaard. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 
// 1. Redistributions of source code must retain the above copyright 
//    notice, this list of conditions and the following disclaimer.  
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.  
// 3. Neither the name of the project nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission. 
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF 
// SUCH DAMAGE.
// 
import "math";

const FLT_RADIX    = 2
const FLT_ROUNDS   = 1
const FLT_DIG      = 6
const FLT_EPSILON  = 1.192092896e-07
const FLT_MANT_DIG = 24
const FLT_MAX      = 3.402823466e+38
const FLT_MAX_EXP  = 38
const FLT_MIN      = 1.175494351e-38
const FLT_MIN_EXP  = (-37)
const DBL_DIG      =  15
const DBL_EPSILON  =  2.2204460492503131e-016
const DBL_MANT_DIG =  53
const DBL_MAX      =  1.7976931348623158e+308
const DBL_MAX_EXP  =  308
const DBL_MIN      =  2.2250738585072014e-308
const DBL_MIN_EXP  =  (-307)
const HUGE_VAL = DBL_MAX

func isdigit(c byte) bool { return c >= '0' && c <= '9' }

func Strtod(str string) float64 {
    // Special-case NaNs
    if str[0] == 'N' { return math.NaN() }
    var number, p10 float64;
    var exponent, n, num_digits, num_decimals int;
    var negative bool;
    slen := len(str);
    
    // char *p = (char *) str;
    p := 0;

    // Skip leading whitespace
    // while (isspace(*p)) p++;
    for ; str[p] == ' '; p++ { }

    // Handle optional sign
    negative = false;
    /*  switch (str[p])
    {             
    case '-': negative = 1; // Fall through to increment position
    case '+': p++;
    } */
    if str[p] == '-' {
        negative = true;
        p++;
    }
    if str[p] == '+' { p++ }
  
    // Process string of digits
    for ; p < slen && isdigit(str[p]) ; {
        number = number * 10. + float64(str[p] - '0');
        p++;
        num_digits++;
    }

    // Process decimal part
    if p < slen && str[p] == '.' {
        p++;
        for ; p < slen && isdigit(str[p]); {
            number = number * 10. + float64(str[p] - '0');
            p++;
            num_digits++;
            num_decimals++;
        }
        exponent -= num_decimals;
    }

    if num_digits == 0 {
        // errno = ERANGE;
        return 0.0;
    }

    // Correct for sign
    if negative { number = -number }

    // Process an exponent string
    if p < slen && (str[p] == 'e' || str[p] == 'E') {
        // Handle optional sign
        negative = false;
        p++;
        if str[p] == '-' {
            negative = true;
            p++;
        }
        if str[p] == '+' { p++ }

        // Process string of digits
        n = 0;
        for ; p < slen && isdigit(str[p]) ; {   
          n = n * 10 + int(str[p] - '0');
          p++;
        }

        if negative {
            exponent -= n
        } else {
            exponent += n
        }
    }

/*    if exponent < DBL_MIN_EXP  || exponent > DBL_MAX_EXP {
        // errno = ERANGE;
        return HUGE_VAL;
    }
*/  
    // Scale the result
    p10 = 10.;
    n = exponent;
    if n < 0 { n = -n }
    for ; n != 0 ; {
        if n & 1 == 1 {
            if exponent < 0 {
                number /= p10
            } else {
                number *= p10;
            }
        }
        n >>= 1;
        p10 *= p10;
    }

  // if number == HUGE_VAL) errno = ERANGE;
  return number;
}
