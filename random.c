/* File: random.c Version 4.2, April 28, 2024 */

#include <stdio.h>
#include <limits.h>  // For ULONG_MAX
#include <time.h>    // For time() for the seed
#include <math.h>    // For log and log2
#include <stdbool.h> // For boolean

/****************************************************************************************
 * Seed generator (randomizer) for all pseudo-random number generators
 * Uses time to generate the seed
 ****************************************************************************************/
unsigned long int state;
void randomize(void) { state = time(NULL); }

/****************************************************************************************
 * unsigned long int pseudo-random number generator with uniform distribution
 * uses xorshift64 algorithm
 * Xorshift random number generators, also called shift-register generators, were invented by George Marsaglia.
 * They generate the next number in their sequence by repeatedly taking the exclusive or of a number with a bit-shifted version of itself.
 * This makes execution extremely efficient on modern computer architectures, but it does not benefit efficiency in a hardware implementation.
 * Like all LFSRs, the parameters have to be chosen very carefully in order to achieve a long period.
 * For execution in software, xorshift generators are among the fastest PRNGs, requiring very small code and state.
 * However, they do not pass every statistical test without further refinement.
 * This weakness can be amended by combining them with a non-linear function, as described in the original paper.
 * The algorithm taken here is A_1(12, 25, 27)·M_{32} which is the best xorshift64 suggested in
 * An experimental exploration of Marsaglia’s xorshift generators, scrambled
 * Sebastiano Vigna, Universit`a degli Studi di Milano, Italy
 ****************************************************************************************/
unsigned long int xorshift64(void) { state ^= state >> 12;
                                     state ^= state << 25;
                                     state ^= state >> 27;
                                     return state * 0x2545F4914F6CDD1DUL;
}

/****************************************************************************************
 * cointoss  gives the result of a coin toss (hence the name) with a theoretical equal probability of getting both sides
 * rantrit   gives the result of a coin of three sides toss with a theoretical equal probability of getting any of the three sides
 * ranfour   gives the result of a dice of four sides toss with a theoretical equal probability of getting any of the four sides
 * uchar_ran gives a random unsigned char
 * uint_ran  gives a random unsigned integer
 * ulong_limited_ran gives a random unsigned int in the interval [0, lim) = [0, lim-1]
 *                   The case lim = 0 has been reset to give a full ulong_ran
 *                   When lim = 1 it returns a non-random 0
  ****************************************************************************************/
bool              cointoss(void)  { unsigned char r = (xorshift64() & 1UL); return r; }
unsigned char     rantrit(void)   { unsigned long r = xorshift64(); return ((r > 6148914691236517205UL) ? ((r > 12297829382473034410UL) ? 2 : 1) : 0); }
unsigned char     ranfour(void)   { unsigned long r = xorshift64(); return ((r < 9223372036854775808UL) ? ((r < 4611686018427387904UL) ? 0 : 1) : ((r < 13835058055282163712UL) ? 2 : 3)); }
unsigned char     uchar_ran(void) { unsigned char r = (xorshift64() & 0xFFUL); return r; }
unsigned int      uint_ran(void)  { unsigned int  r = (xorshift64() & 0xFFFFFFFFUL); return r; }
unsigned long int ulong_limited_ran (const unsigned long int lim) { unsigned long int r;
                  /*switch (lim) {
                      case 0: return xorshift64(); break;
                      case 1: return 0; break;
                      case 2: return cointoss(); break;
                      case 3: return rantrit();  break;
                      case 4: return ranfour();  break;
                     default: r = (xorshift64() % lim); return r; break; }*/ // unnecessary steps
                     if (lim == 0) return xorshift64();  // No limit case
                        return xorshift64() % lim;
}

/*
 generates a random double in the interval [min,max] with a step
*/
double discrete_random(int power, double min,double step) {
    unsigned long u = ulong_limited_ran((1<<power)); 
    return min + u * step;
}
/****************************************************************************************
 * 2 sided coin toss with asymmetric probability p 1-p aka bernoulli_trial
 *
 * The double pseudo-random number generator between zero and one with uniform
 * distribution based on xorshift64 function is mal rotllo because the mantissa
 * of a double has 52 bits (one implicit). So, 12 bits of the integer random number
 * are lost.
 * double double_uniform(void){ return ldexp(((double) xorshift64()), -64); }
 * It is better to check whether we got the price in the lottery without transforming
 * the output of xorshift64().
 * NOTE: the real number in [0,1) with uniform distribution associated to xorshift64()
 * is xorshift64() / 2^{64} = xorshift64() * 2^{-64}.
 *****************************************************************************************/
 bool bernoulli_trial (const double prob) { return (xorshift64() < ldexp(prob, 64)); }

/****************************************************************************************
    This function computes the value of log(1+x) in a way that is accurate for small x
*****************************************************************************************/
double gsl_log1p (const double x) { volatile double y, z;
       y = 1 + x;
       z = y - 1;
       return log(y) - (z-x)/y ;  /* cancels errors with IEEE arithmetic */
}

/*************************************************************************************************************************************
 * double pseudo-random number generator with exponential distribution
 * GSL dixit:  double u = gsl_rng_uniform (r); return -mu * log1p (-u);
 * EXPLANATION:
 * Exponential distribution with mean \mu:        f(x,\mu) = exp(-x/\mu) / \mu for x >= 0 & 0 for x < 0
 *    The cumulative exponential distribution is: F(x,\mu) = 1 - exp(-x/\mu) for x >= 0 & 0 for x < 0
 * Consecuently a sample from the exponential distribution with mean \mu can be obtained as:
 *      F^{-1}(u,\mu) where u ~ U[0,1) <==> u = 1 - exp(-x/mu) <==>  x = -mu * log_e(1-u)
 * **************************************************************************************************************************
 * Since we need to use the exponential distribution with mean \mu = 1 and we need to compute (integer) number of occurrences
 * we must compute         floor(-log (1 - x / 2^{64}))               where x = xorshift64()
 * NOTA LLUIS: De fet hauria de ser Round_to_nearest_integer(-log (1 - x / 2^{64}))
 * THE SAME COMMENT AS BEFORE APPLIES: It is not good to take double_uniform = xorshift64() / 2^{64} as a real number in [0,1) with uniform distribution.
 *
 * A MORE FRIENDLY FORMULA for floor(-log (1 - x / 2^{64})): We define Ψ[y] := floor(-log((1 + y) / 2^{64})
 * Then, since x = xorshift64() \in [0, ULONG_MAX] \cap \Z ===> 2^{64} - x \in [1, 2^{64}] \cap \Z,
 *     Ψ[ULONG_MAX - x] = floor(-log((1 + ULONG_MAX - x) / 2^{64}) = floor(-log((2^{64} - x) / 2^{64}) = floor(-log (1 - x / 2^{64})) = floor(64*log(2) - log(2^{64} - x))
 *
 * FACTS:
 * *   0: -log(1 - x / 2^{64}) is increasing for x \in [0, ULONG_MAX]
 * *   I: -log(1 - x / 2^{64}) >= -log (1 - x / 2^{64})\evalat{x=0} = -log (1) = 0
 * *  II: -log (1 - 11660566172440666342UL / 2^{64}) > 1 && -log (1 - 11660566172440666341UL / 2^{64}) < 1
 *        Consequently x < 11660566172440666342UL <===> 0 <= -log (1 - x / 2^{64}) < 1
 * * III: -log (1 - 15950248739700762817UL / 2^{64}) > 2 && -log (1 - 15950248739700762816UL / 2^{64}) < 2
 *        Consequently x < 15950248739700762817UL <===> 0 <= -log (1 - x / 2^{64}) < 2 &
 * *  IV: floor(-log (1 - (x = (ULONG_MAX - [k=0:27])) / 2^{64}))   = Ψ[k=0:27]   = expWP1[k]
 *        floor(-log (1 - (x = (ULONG_MAX - [k=28:77])) / 2^{64}))  = Ψ[k=28:77]  = 40
 *        floor(-log (1 - (x = (ULONG_MAX - [k=78:212])) / 2^{64})) = Ψ[k=78:212] = 39
 *        Ψ[213] = 38
 * *   V: Assume that k = ULONG_MAX - x >= 213
 *        We have 1 + k = 1 + ULONG_MAX - x = 2^{64} - x >= 214 and
 *                floor(-log (1 - x / 2^{64})) = Ψ[k] = floor(64*log(2) - log(1 + k))
 *        On the other hand, 64*log(2) = 44.361419555...., and
 *                           x >= 15950248739700762817UL ==> log(1 + k) = log(2^{64} - x) <= log(2496495334008788799) = 42.361419555....
 *        So, there are no cancellations in the computation of 64*log(2) - log(1+k)
 *************************************************************************************************************************************/
#define SIXTYFOUR_TIMES_LOGTWO 44.36141955583649980270285577332330035683
unsigned int exponential_with_parameter_one(void) { unsigned long int x = xorshift64();
             if (x < 15950248739700762817UL) return ((x < 11660566172440666342UL) ? 0 : 1);
                                       /*  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27 */
             unsigned char expWP1[28] = { 44, 43, 43, 42, 42, 42, 42, 42, 42, 42, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41 };
             unsigned long int k = ULONG_MAX - x; if (k < 213) return ((k < 28) ? expWP1[k] : ((k < 78) ? 40 : 39));
             return ((unsigned int) (SIXTYFOUR_TIMES_LOGTWO - log1p(((double) k))));
}

/*******************************************************************************************************************************************************
 * The exponential distribution with mean \mu truncated at the interval [0,T] is f_{[0,T]}(x,\mu) = \frac{(exp(-x/\mu) / \mu) * \Xi([0,T])}{C}
 * The cumulative exponential distribution with mean \mu truncated at the interval [0,T] is
 *     F_{[0,T]}(x,\mu) = \int_{0}^{x} f_{[0,T]}(s,\mu) ds = \int_{0}^{x} \frac{(exp(-s/\mu) / \mu) * \Xi([0,T])}{C} ds
 *                      = \int_{0}^{x} exp(-s/\mu) / (C * \mu) ds = (\mu - \mu * exp(-x/\mu)) / (C*\mu),
 * and the value of C can be determined with F_{[0,T]}(T,\mu) = 1 which gives
 *     C = 1 - exp(-T/mu) = F(T,\mu)       and        T > 0
 * When the mean is 1 (that is, \mu = 1) this amounts F_{[0,T]}(x,1) = (1 - exp(-x)) / (1 - exp(-T)).
 * Consecuently, as above, a sample from the exponential distribution with mean 1 truncated at the interval [0,T] can be obtained as:
 *      F^{-1}_{[0,T]}(u,1) where u ~ U[0,1) <==> u = (1 - exp(-x)) / (1 - exp(-T)) <==>
 *      u * (1 - exp(-T)) = 1 - exp(-x) <==> exp(-x) = 1 - u * (1 - exp(-T)) <==> x = -log_e(1 - u * (1 - exp(-T)))
 * Since we need to compute (integer) number of occurrences we must compute
 *                   floor(-log (1 - z * (1 - exp(-T)) / 2^{64}))               where z = xorshift64() \in [0,ULONG_MAX]
 * NOTA LLUIS: De fet hauria de ser Round_to_nearest_integer(-log (1 - z * (1 - exp(-T)) / 2^{64}))
 * *******************************************************************************************************************************************
 * VERY IMPORTANT OBSERVATION:
 * -log(1 - z * (1 - exp(-T)) / 2^{64}) is increasing as a function of z.
 * Thus, -log (1 - z * (1 - exp(-T)) / 2^{64}) < -log (1 - 2^{64} * (1 - exp(-T)) / 2^{64}) = -log (1 - (1 - exp(-T))) = -log (exp(-T)) = T
 * In particular floor(-log (1 - z * (1 - exp(-T)) / 2^{64})) = 0 whenever T = 1;
 *
 * VERY IMPORTANT OBSERVATION II:
 * -log (1 - 11660566172440666341 * (1 - exp(-T)) / 2^{64})) < -log (1 - 11660566172440666341 / 2^{64})) < 1
 * Thus, for z <= 11660566172440666341 and every T >= 1, floor(-log (1 - z * (1 - exp(-T)) / 2^{64})) = 0.
 * *******************************************************************************************************************************************
 * for T:2 thru 20 do ( for s:1 thru T-1 do ( aa: solve(-log (1 - z * (1 - exp(-T)) / 2^(64)) = s,z), print(T, s, bfloat(aa), floor(aa))) )$
 * *******************************************************************************************************************************************/
/*unsigned short int truncated_exponential_with_parameter_one (unsigned short trunc){ unsigned long int z = xorshift64();
         switch (trunc) {
            case 10: { unsigned long int zboundaries[9] = { 11661095585361195362UL, 15950972912750644994UL, 17529130586644716243UL,
                                                            18109702349797290482UL, 18323282765585778462UL, 18401854609591211717UL,
                                                            18430759575655740243UL, 18441393118418638501UL, 18445304980187926148UL };
                       for (register int t=0; t < 9 ; t++) { if (z < zboundaries[t]) return (t+1); }
                       return 10;
                     } ; break;

             case 9: { unsigned long int zboundaries[8] = { 11662005378239652701UL, 15952217399724810975UL, 17530498200695774902UL,
                                                            18111115259768589493UL, 18324712338994902839UL, 18403290313136531285UL,
                                                            18432197534352137602UL, 18442831906738754116UL };
                       for (register int t=0; t < 8 ; t++) { if (z < zboundaries[t]) return (t+1); }
                       return 9;
                     } ; break;
             case 8: { unsigned long int zboundaries[7] = { 11664479169275905745UL, 15955601247625055715UL, 17534216839806579139UL, 18114957061682843469UL,
                                                            18328599449972463014UL, 18407194092386980546UL, 18436107445517502593UL };
                       for (register int t=0; t < 7 ; t++) { if (z < zboundaries[t]) return (t+1); }
                       return 8;
                     } ; break;
             case 7: { unsigned long int zboundaries[6] = { 11671208937386888134UL, 15964806759067919697UL, 17544333126322859647UL,
                                                            18125408403624165328UL, 18339174051916310545UL, 18417814039151676011UL };
                       for (register int t=0; t < 6 ; t++) { if (z < zboundaries[t]) return (t+1); }
                       return 7;
                     } ; break;
             case 6: { unsigned long int zboundaries[5] = { 11689541649247971805UL, 15989883698723615902UL, 17571891128730770820UL,
                                                            18153879138010872624UL, 18367980561633316663UL };
                       for (register int t=0; t < 5 ; t++) { if (z < zboundaries[t]) return (t+1); }
                       return 6;
                     } ; break;
             case 5: { unsigned long int zboundaries[4] = { 11739667429366628617UL, 16058449722820606775UL, 17647240939477576282UL, 18231724564399438168UL };
                       for (register int t=0; t < 4 ; t++) { if (z < zboundaries[t]) return (t+1); }
                       return 5;
                     } ; break;
             case 4: return ((z < 11878121557565525290UL) ? 1 : ((z < 16247838278329192791UL) ? 2 : ((z < 17855367223641238263UL) ? 3 : 4))); break;
             case 3: return ((z < 12271529658528073609UL) ? 1 : ((z < 16785973131626161745UL) ? 2 : 3)); break;
             case 2: return ((z < 13485650502877570763UL) ? 1 : 2); break;
             case 0: case 1: return 1; break;
            default: if (z < 11660566172440666342UL) return 1; else  return (((unsigned short int) - log1p (-ldexp(z * (1.0 - exp(-trunc)), -64))) + 1);
                     break;
         }
}*/
/****************************************************************************************
 * unsigned integer pseudo-random number generator with geometric distribution
 *
 * double u = gsl_rng_uniform_pos (r); if (p == 1) k = 1; else k = log (u) / log (1 - p) + 1; return k;
 * p must be positive. If p <= 0 the return value is 0 (to be justified later), which is theoretically impossible.
 *
 * SAME COMMENT AS BEFORE: It is not good to take double_uniform = gsl_rng_uniform_pos (r) = xorshift64() / 2^{64}
 * as a real number in [0,1) with uniform distribution.
 *
 * Since u \in [0, 1) ===> log(u) is increasing with u and, hence, log(u) < log(1) = 0.
 * floor(log(u) / log(1-p)) = n <===> n <= log (u) / log (1 - p) < n+1
 * JUSTIFICATION OF THE RETURN VALUE WHEN p <= 0: 'if (p <= 0) return 0;'
 * If p < 0 ===> log (1 - p) > log(1) = 0. Hence, floor(log(u) / log(1-p)) = n <===>
 *    n * log (1 - p) <= log (u) < (n+1) * log (1 - p) ===> n <= -1 ===> 1 + floor(log(u) / log(1-p)) <= 0
 * If p = 0 ===> log (u) / log (1 - p) = -\infty < (n+1) * log (1 - p) for every n ===> 1 + floor(log(u) / log(1-p)) = -\infty < 0
 *
 * ANOTHER LIMIT CASE p = 1: 'if (p >= 1) return 1;'
 * Assume that u \in (0,1) and take p \in (1-u,1) so that 1-p < u. Clearly
 *        log(1-p) < log(u) < 0  and hence   0 <= log (u) / log (1 - p) < 1
 * Consequently, 1 + log(u) / log(1-p) \in [1, 2). Moreover,
 *        lim_{p\to 1} (1 + log(u) / log(1-p)) \in [1, 2)
 *        lim_{u\to 0 ===> p\to 1} (1 + log(u) / log(1-p)) \in [1, 2)
 * Hence, lim_{p\to 1} 1 + floor(log(u) / log(1-p)) = 1 for every u \in [0, 1)
 *
 * From now on we assume p, u \in (0, 1): Recall that log\evalat{(0,1)} is strictly increasing and negative.
 * Then, floor(log(u) / log(1-p)) = n <===> n <= log (u) / log (1 - p) < n+1 <===>
 *        n*log (1 - p) >= log (u) > (n+1) * log (1 - p) <===> log ((1 - p)^n) >= log (u) > log ((1 - p)^{n+1})
 *        <===> (1-p)^n >= u > (1-p)^{n+1}
 * Write u = x / 2^{64} with x = xorshift64() \in [1, ULONG_MAX] \cap \Z. Then,
 *       floor(log(u) / log(1-p)) = n <===> 2^{64} * (1-p)^n >= x > 2^{64} * (1-p)^{n+1}
 * Observe that x = 0 < 2^{64} * (1-p)^n for every n ===> 1 + floor(log(x=0 / 2^{64}) / log(1-p)) = \infty
 *
 * Condition for return 1 <===> floor(log(u) / log(1-p)) = 0 <===> 2^{64} >= x > 2^{64} * (1-p)
 *           Observe that since 2^{64} * (1-p)^0 = 2^{64} > ULONG_MAX >= x the above condition is equivalent to x > 2^{64} * (1-p)
 * Condition for return 2 <===> floor(log(u) / log(1-p)) = 1 <===> 2^{64} * (1-p) >= x > 2^{64} * (1-p)^2
 *******************************************************************************************************************/
unsigned int geometric (const double p) { unsigned int k,x;
/* First we deal with special and wrong parameters, and limit cases */
             if (p >= 1) return 1;
             if (p <= 0) return 0;
             if ((x = xorshift64()) == 0) return UINT_MAX;

             double pp = 1-p, bd = ldexp(pp, 64), bd2 = bd * pp;
             if ( x > bd2 * pp) return ((x > bd) ? 1 : ((x > bd2) ? 2 : 3));

             double a = - gsl_log1p(-p);
             k = ((a + SIXTYFOUR_TIMES_LOGTWO) - log(xorshift64())) / a;
             return k;
}