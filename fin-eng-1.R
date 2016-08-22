# @rob4lderman
# July 2016
#
# Functions for building binomial models, asset price lattices,
# bond payment schedules, mortgage/MBS payment schedules, CDS payment schedules.
#
#   calibrateBinomialModel 
#   buildAssetPriceLattice
#   buildShortRateLattice
#
#   buildBDTShortRateLattice 
#   calibrateBDTObjectiveFn 
#   calibrateBDTShortRateLattice 
#
#   buildElementarySecurityPriceLattice 
#   buildFuturesPriceLattice 
#
#   buildInterestRateSwapCouponPaymentLattice 
#   buildInterestRateSwapPriceLattice 
#
#   buildOptionPriceLattice 
#   shouldExerciseOption 
#   shouldExerciseOptionEarly 
#
#   buildZCBPriceLattice 
#   buildDefaultableZCBPriceLattice 
#   buildBondPriceLattice 
#
#   buildBinomialPriceLattice 
#
#   buildDefaultableBondSchedule 
#   computeDefaultableBondPricing 
#
#   buildCDSSchedule 
#
#   computeMonthlyMortgagePayment 
#   computeMortgagePresentValue 
#   computeMortgagePrincipalPresentValue 
#   buildMortgageSchedule 
#
#   buildCPRSchedule 
#   buildMBSPassThruSchedule 
#   buildMBSTrancheSchedule 
#

#
# Calibrate a binomial model from the given Black-Scholes parameters.
#
# @param T: time to maturity (in years)
# @param n: number of periods in the binomial model
# @param r: annualized continuously compounded market interest rate
# @param c: annualized (continuously compounded?) dividend yield
# @param sigma: annualized volatility of the underlying asset
#
# @return a list object: list(q=q, u=u, d=d, R=R, n=n, r=r, sigma=sigma, T=T, c=c)
#           q: risk-neutral probability of an upmove (1-q is the risk-neutral 
#              probabilty for a downmove)
#           u: magnitude/return of an upmove
#           d: magnitude/return of a downmove
#           R: periodically compounded return (note: R = 1 + periodic-rate-of-return )
#           n: number of periods in the model
#         Note: q, u, d, R, n are the the only fields relevant to the binomial model.
#         The rest are included for completeness.
#
calibrateBinomialModel <- function(T,
                                   n,
                                   r,
                                   c,
                                   sigma) {

    R = exp(r * T/n)        # convert annualized continuously-compounded rate r into periodic-compounded return R
    C = exp(-c * T/n)       # convert annualized continuously-compounded yeild c into a periodic-compounded return C

    u = exp(sigma * sqrt(T/n))  # convert annualized volatility into periodic volatility
    d = 1 / u

    q = (R * C - d) / (u - d)   # risk/volatility-neutral probability of an upmove

    list(q=q, u=u, d=d, R=R, n=n, r=r, sigma=sigma, T=T, c=c)
}


#
# Construct an asset price lattice where in every period the asset either 
# makes an upmove by a factor of u or a downmove by a factor of d.
# 
# @param S0: asset spot price at t=0
# @param n:  number of periods in the binomial model
# @param u:  upmove factor
# @param d:  downmove factor (default = 1/u)
# 
# @return asset.price.lattice: n+1 x n+1 matrix for t=0..n
#         only the upper half of the matrix is filled in (the rest = NA)
#      
#         [ S_0   S_0*u^1*d^0     S_0*u^2*d^0   S_0*u^3*d^0 ... ]
#         [  .    S_0*u^0*d^1     S_0*u^1*d^1   S_0*u^2*d^1 ... ]
#         [  .         .          S_0*u^0*d^2   S_0*u^1*d^2 ... ]
#         [  .         .               .        S_0*u^0*d^3 ... ]
#         [  .         .               .             .      ... ]
#
#           t=0       t=1             t=2           t=3     ...
#
#         Note that the j indices (states of the world) are reversed:
#
#           asset.price.lattice[j,t] corresponds to S_t-1,t-j
#
#           asset.price.lattice[1,1] corresponds to S_0,0
#
#           asset.price.lattice[1,2] corresponds to S_1,1
#           asset.price.lattice[2,2] corresponds to S_1,0
#
buildAssetPriceLattice <- function(S0,
                                   n,
                                   u,
                                   d = 1/u) {     
                                   
    asset.price.lattice <- matrix(NA, nrow=n+1, ncol=n+1)    

    # for each time period 
    for (i in 1:(n+1)) {            # iterating over cols
        t = i - 1                   # t=0..n, i=1..n+1
        for (j in 0:t) {            # iterating over rows
            asset.price.lattice[j+1,i] = S0 * u^(t-j) * d^(j)
        }
    }
    asset.price.lattice
}


#
# Construct a short-rate lattice where in every period the interest rate 
# either makes an upmove by a factor of u or a downmove by a factor of d.
#
# Note: this method maps directly to buildAssetPriceLattice.
# 
# @param r0: short rate at t=0
# @param n:  number of periods in the binomial model
# @param u:  upmove factor
# @param d:  downmove factor (default = 1/u)
# 
# @return short.rate.lattice: n+1 x n+1 matrix for t=0..n
#         only the upper half of the matrix is filled in (the rest = NA)
#      
#         [ r0   r0*u^1*d^0     r0*u^2*d^0   r0*u^3*d^0 ... ]
#         [  .   r0*u^0*d^1     r0*u^1*d^1   r0*u^2*d^1 ... ]
#         [  .         .        r0*u^0*d^2   r0*u^1*d^2 ... ]
#         [  .         .            .        r0*u^0*d^3 ... ]
#         [  .         .            .            .      ... ]
#
#           t=0       t=1         t=2           t=3     ...
#
#         Note that the j indices (states of the world) are reversed:
#
#           short.rate.lattice[j,t] corresponds to r_t-1,t-j
#
#           short.rate.lattice[1,1] corresponds to r_0,0
#
#           short.rate.lattice[1,2] corresponds to r_1,1
#           short.rate.lattice[2,2] corresponds to r_1,0
#
buildShortRateLattice <- function(r0,
                                  n,
                                  u,
                                  d) {
    buildAssetPriceLattice(r0,n,u,d)
}


#
# Build a BDT short-rate lattice using the following formula:
#
#     r_i,j = a_i * exp(b_i * j)
#
# @param a: vector of a_i values
# @param b: vector of b_i values
# @param n: the number of time periods in the lattice (default: length(a) - 1)
#
# @return a BDT short rate lattice: n+1 x n+1 matrix, for time periods 0 thru n
#   
#         Note that the j indices (states of the world) are reversed:
#
#           bdt.short.rate.lattice[j,t] corresponds to r_t-1,t-j
#
#           bdt.short.rate.lattice[1,1] corresponds to r_0,0
#
#           bdt.short.rate.lattice[1,2] corresponds to r_1,1
#           bdt.short.rate.lattice[2,2] corresponds to r_1,0
#
buildBDTShortRateLattice <- function(a,
                                     b,
                                     n = length(a)-1) {
                                   
    bdt.short.rate.lattice <- matrix(NA, nrow=n+1, ncol=n+1)    

    # for each time period 
    for (i in 0:n) {                # iterating over cols (time periods)
        for (j in 0:i) {            # iterating over rows (states of the world)
            bdt.short.rate.lattice[j+1,i+1] = a[i+1] * exp(b[i+1] * (i-j))
        }
    }

    bdt.short.rate.lattice
}



#
# Optimization objective function used with calibrateBDTShortRateLattice
#
# 1. Builds a BDT short-rate lattice using the given a,b,n.
# 2. Builds an elem price lattice using the BDT short-rate lattice, n+1, and q
# 3. Computes ZCB prices from the elem price lattice
# 4. Computes spot rates from the ZCB prices
# 5. Compares model spot rates (4) with market spot rates (s)
#
# @param a: vector of a_i values (which are being tuned)
# @param b: vector of b_i values
# @param s: vector of market spot rates
# @param q: risk-neutral probabilities
# @param n: number of time periods in the BDT short rate lattice
#
# @return sum of squared diffs between market spot rates and model-computed spot rates
#
calibrateBDTObjectiveFn <- function(a, b, s, q, n) {

    bdt.short.rate.lattice <- buildBDTShortRateLattice(a,b,n)

    elem.price.lattice <- buildElementarySecurityPriceLattice(short.rate.lattice=bdt.short.rate.lattice,
                                                              n=n+1,    # one more time period 
                                                              q=q)

    bdt.zcb.prices <- colSums(elem.price.lattice[,2:(n+2)], na.rm=T)

    # zcb.price = 1/(1+s_i)^i
    # 1/zcb.price = (1 + s_i)^i
    # (1/zcb.price)^(1/i) - 1 = s_i
    bdt.spot.rates <- (1/bdt.zcb.prices)^(1/1:(n+1)) - 1

    sq.diffs <- (s*100 - bdt.spot.rates*100)^2
    sum(sq.diffs)       # to minimize
}


#
# Calibrate a BDT short-rate lattice using the given market spot rates.
#
# The a_i values are calibrated using R's optim() function.
#
# @param a: vector of initial a_i values (to be calibrated)
# @param b: vector of b_i values
# @param s: vector of market spot rates
# @param q: risk-neutral probabilities. default: 0.5
# @param n: number of time periods in the calibrated BDT short rate lattice. default: length(a)-1
#
# @return list: [ $par: calibrated a_i values (see ?optim),
#                 $value: final value of calibrateBDTObjectiveFn (see ?optim),
#                 $counts: see ?optim,
#                 $convergence: see ?optim,
#                 $message: see ?optim,
#                 $bdt: calibrated BDT short-rate lattice 
#               ]
#
calibrateBDTShortRateLattice <- function( a,
                                          b,
                                          s,
                                          q = 0.5,
                                          n = length(a)-1,
                                          optim.rounds = 1) {


    for (i in 1:optim.rounds) {
        o <- optim(par=a,
                   fn=calibrateBDTObjectiveFn,
                   b=b,
                   s=s,
                   q=q,
                   n=n,
                   method = "L-BFGS-B" )
                   # control = list(trace=1))
        a <- o$par
    }

    o$bdt <- buildBDTShortRateLattice(a, b, n)
    o
}


#
# An "elementary security" pays $1 at any given node and $0 everywhere else.
#
# Each node in the elementary security price lattice gives the t=0 price of
# an elementary security, where the elementary security pays $1 at the given 
# node and $0 at every other node.
# 
#
#             1       Pe_k,s-1        1        Pe_k,s
# Pe_k+1,s = --- * --------------- + --- * ----------------, for 0 < s < k+1
#             2    (1 + r_k,s-1)      2      (1 + r_k,s)
# 
#             1       Pe_k,0
# Pe_k+1,0 = --- * -------------      <--- the 1/2 comes from risk-neutral probabilties (q=0.5)
#             2     (1 + r_k,0)       <--- discounted
# 
#               1       Pe_k,k
# Pe_k+1,k+1 = --- * ------------      
#               2     (1 + r_k,k)     <-- discounted
# 
# 
# Pe_0,0 = 1              <--- price at t=0 of $1 paid at t=0
#
# @param short.rate.lattice: interest rates for each node
# @param n: number of time periods in the price lattice
# @param q: risk-neutral probability (default: 0.5)
#
# @return elem.price.lattice: n+1 x n+1 matrix
#
buildElementarySecurityPriceLattice <- function(short.rate.lattice,
                                                n=ncol(short.rate.lattice)-1,
                                                q=0.5) {

    elem.price.lattice <- matrix(NA, nrow=n+1, ncol=n+1)    

    elem.price.lattice[1,1] = 1

    # for each time period 
    for (t in 2:(n+1)) {            # iterating over cols (time periods)
        for (j in 1:t) {            # iterating over rows (states of the world)

            # upmove <- bin.price.lattice[j,t+1]
            # downmove <- bin.price.lattice[j+1,t+1]

            if (j == 1) {
                # top branch of lattice (upmove only)
                r <- short.rate.lattice[j,t-1]
                elem.price.lattice[j,t] <- q * elem.price.lattice[j,t-1] / (1 + r)

            } else if (j == t) {
                # bottom branch of lattice (downmove only)
                r <- short.rate.lattice[j-1,t-1]
                elem.price.lattice[j,t] <- (1-q) * elem.price.lattice[j-1,t-1] / (1 + r)

            } else {
                # middle of lattice (upmove and downmove)
                r1 <- short.rate.lattice[j-1,t-1]   # downmove
                r2 <- short.rate.lattice[j,t-1]     # upmove

                elem.price.lattice[j,t] <- (1-q) * elem.price.lattice[j-1,t-1] / (1 + r1) +
                                           q * elem.price.lattice[j,t-1]/ (1 + r2)
            }
        }
    }

    elem.price.lattice 
}


#
#
# A futures price lattice is a binomial price lattice where the 
# final prices (t=n) are the underlying asset's final prices, and
# there is NO DISCOUNTING of the expected values for earlier nodes
# in the lattice.
#
#   S_t = E_Q_t[ S_t+1 ]            <--- note: no discounting
# 
# @param asset.price.lattice: as built by buildAssetPriceLattice()
# @param n:  number of periods in the binomial model.
#            must be <= periods in asset.price.lattice
# @param q:  risk-neutral probability q
# 
# @return futures.price.lattice: n+1 x n+1 matrix for t=0..n
#         only the upper half of the matrix is filled in (the rest = NA)
#
buildFuturesPriceLattice <- function(asset.price.lattice,
                                     n,
                                     final.price.vec = asset.price.lattice[,n+1],
                                     q=0.5) {

    buildBinomialPriceLattice( final.price.vec = final.price.vec,
                               n = n,
                               q = q )
}


#
# The "coupon payments" (i.e intermediate cash flows) for an interest rate swap 
# are simply the difference between the short rate at each node and the "strike rate", K.
#
# @param short.rate.lattice: interest rates at each node
# @param n: number of time periods in the lattice
# @param forward.start: time period at which the swap "activates" (cash flows 
#                       technically begin at t=forward.start+1, but we record
#                       one node earlier and discount)
# @param K: the "strike rate"
#
# @return the coupon.payment.lattice: n+1 x n+1 matrix of coupon payments for each node.
# 
buildInterestRateSwapCouponPaymentLattice <- function( short.rate.lattice,
                                                       n=ncol(short.rate.lattice)-1,    
                                                       forward.start=0,
                                                       K) {

    # Coupon payment for time period t is technically made in time period t+1.
    # However we record it in the coupon payment lattice at time period t.
    # When pricing the swap for time period t, we discount the coupon payment
    # by one period (since it's not paid until t+1)
    coupon.payment.lattice <- matrix(NA, nrow=n+1, ncol=n+1)

    coupon.payment.lattice[1:(n+1),1:(n+1)] <- short.rate.lattice[1:(n+1),1:(n+1)] - K

    # Handle "forward-start" swaps: swaps that don't begin paying coupons
    # until t=forward.start (technically forward.start+1, but we record one
    # node earlier and discount).  Set all previous coupons in the lattice to 0.
    if (forward.start > 0) {
        for (t in 1:forward.start) {
            for (j in 1:t) {
                coupon.payment.lattice[j,t] = 0
            }
        }
    }

    coupon.payment.lattice
}



#
# An interest rate swap price lattice is a binomial price lattice where the final
# prices (t=n) are equal to the "coupon payments" for those nodes.
#
# Note that interest-rate swap "coupons" are paid in arrears, therefore all coupon
# payments are discounted, just like the expected values.
# 
# @param short.rate.lattice: the interest rate at each node
# @param n: number of time periods in the lattice
# @param forward.start: time period at which the interest rate swap "activates" (cash flows begin)
# @param q: risk-neutral probabilities
# @param K: "strike rate"
#
# @return interest rate swap price lattice: n+1 x n+1 matrix of the swap price at each node.
#
buildInterestRateSwapPriceLattice <- function( short.rate.lattice,
                                               n=ncol(short.rate.lattice)-1,    
                                               forward.start=0,
                                               q = 0.5,
                                               K) {

    coupon.payment.lattice <- buildInterestRateSwapCouponPaymentLattice(short.rate.lattice,
                                                                        n=n,
                                                                        forward.start=forward.start,
                                                                        K=K)
    
    buildBinomialPriceLattice( short.rate.lattice = short.rate.lattice,
                               n = n,
                               q = q,
                               coupon.payment.lattice = coupon.payment.lattice,
                               discount.coupons = T )
}


#
# Construct an option price lattice based on the given asset.price.lattice.
#
# The t=n prices in the asset.price.lattice are used to determine the 
# option payoff at maturity.  
#
# The option price lattice is filled in by working backwards thru the lattice, 
# computing prices at each node as the discounted expected value (under Q probability 
# space) of the option in the next time step.
# 
#           1
#   S_t = ----- * E_Q_t[ S_t+1 ]
#          r_t
#
#           1
#       = ----- * (q * upmove-payoff + (1-q) * downmove-payoff)
#          r_t 
#       
# If an American option, the intrinsic value is computed at each node and
# compared to the discounted expected payoff at that node.  The option price
# (value) at the node is the max(intrinsic-value, discounted-expected-payoff).
# If the intrinsic value is greater than the discounted expected payoff, 
# this indicates a node where it's optimal to exercise the option early.
#
# @param asset.price.lattice: as built by buildAssetPriceLattice()
# @param n:  number of periods in the binomial model.
#            must be <= periods in asset.price.lattice
# @param K:  option strike price
# @param R:  1-period return (optional - if not specified, short.rate.lattice 
#            must be specified)
# @param short.rate.lattice: short rates for each time period and state 
#                            (optional - if not specified, R must be specified)
# @param q:  risk-neutral probability q
# @param is.call.option:  T if call option, F if put option
# @param is.amer.option:  T if American option, F if European option
# 
# @return option.price.lattice: n+1 x n+1 matrix for t=0..n
#         only the upper half of the matrix is filled in (the rest = NA)
#
buildOptionPriceLattice <- function(asset.price.lattice,
                                    n,
                                    K,
                                    R,
                                    short.rate.lattice=matrix(R-1,nrow=n,ncol=n),
                                    q=0.5,
                                    is.call.option = T,
                                    is.amer.option = F) {

    # n+1 x n+1 matrix (only the upper half is relevant)
    option.price.lattice <- matrix(NA, nrow=n+1, ncol=n+1)
    
    if (is.call.option) {
        # call option price at maturity (t=n) = max(S_n - K,0)
        option.price.lattice[,n+1] = pmax(asset.price.lattice[1:(n+1),n+1] - K,0)
    } else {
        # put option price at maturity (t=n) = max(K - S_n,0)
        option.price.lattice[,n+1] = pmax(K - asset.price.lattice[1:(n+1),n+1],0)
    }

    for (i in n:1) {            # iterating backwards over cols
        for (j in 1:i) {        # iterating over rows
            op.upmove <- option.price.lattice[j,i+1]
            op.downmove <- option.price.lattice[j+1,i+1]

            r <- short.rate.lattice[j,i]

            op.expected.payoff <- 1 / (1 + r) * (q * op.upmove + (1 - q) * op.downmove)
            op.intrinsic.value <- ifelse(is.call.option,
                                         max(asset.price.lattice[j,i] - K, 0),
                                         max(K - asset.price.lattice[j,i], 0))

            # American options can be exercised early.  
            # Need to check if intrinsic value, max(K-S_t,0), is greater than 
            # expected payoff at maturity.

            option.price.lattice[j,i] = ifelse(is.amer.option,
                                               max(op.expected.payoff, op.intrinsic.value),
                                               op.expected.payoff)
        }
    }

    option.price.lattice
}



#
# Determine whether or not to exercise an option (possibly early),
# based on the diff between its expected payoff and its intrinsic value.
#
# @param asset.price.lattice: as built by buildAssetPriceLattice()
# @param option.price.lattice: as buitl by buildOptionPriceLattice()
# @param K: strike price
# @param is.call.option: T if a call option, F if a put option
#
# @return a boolean matrix indicating, for each node in the binomial model,
#         whether the option's intrinsic value is greater than its expected payoff
#         (and therefore should be exercised early)
#
shouldExerciseOption <- function(asset.price.lattice,
                                 option.price.lattice,
                                      K,
                                      is.call.option = T) {

    option.intrinsic.value <- if (is.call.option) pmax(asset.price.lattice - K,0) else pmax(K - asset.price.lattice,0)

    # only compare with the dimensions of the option.price.lattice
    num.rows <- nrow(option.price.lattice)
    num.cols <- ncol(option.price.lattice)

    option.intrinsic.value[1:num.rows,1:num.cols] == option.price.lattice & option.price.lattice > 0
}
    
#
# Determine whether or not to exercise an option (possibly early),
# based on the diff between its expected payoff and its intrinsic value.
#
# @param asset.price.lattice: as built by buildAssetPriceLattice()
# @param option.price.lattice: as buitl by buildOptionPriceLattice()
# @param K: strike price
# @param is.call.option: T if a call option, F if a put option
#
# @return a boolean matrix indicating, for each node in the binomial model,
#         whether the option's intrinsic value is greater than its expected payoff
#         (and therefore should be exercised early)
#
shouldExerciseOptionEarly <- function(asset.price.lattice,
                                      option.price.lattice,
                                      K,
                                      is.call.option = T) {
    shouldExerciseOption(asset.price.lattice = asset.price.lattice,
                         option.price.lattice = option.price.lattice,
                         K = K,
                         is.call.option = is.call.option)
}


#
#
# Build a zero-coupon bond price lattice.
# 
# @param short.rate.lattice: interest rate at each node
# @param face.value: the face value of the bond
# @param n:  number of periods in the binomial model
# @param q:  risk-neutral probability q
# 
# @return bond.price.lattice: n+1 x n+1 matrix of bond prices at each node.
#
buildZCBPriceLattice <- function(short.rate.lattice,
                                 face.value,
                                 n=ncol(short.rate.lattice)-1,
                                 q=0.5) {

    buildBondPriceLattice( short.rate.lattice = short.rate.lattice,
                           face.value = face.value,
                           coupon.rate = 0,
                           n = n,
                           q = q)
}

#
#
# Build a defaultable zero-coupon bond price lattice.
#
#                    1
# Z'[T]_i,j,0 = ----------- * [ q_u (1 - h_i,j) * Z'[T]_i+1,j+1,0         <-- no default
#                 1 + r_i,j
#                               + q_d (1 - h_i,j) * Z'[T]_i+1,j,0         <-- no default
# 
#                               + q_u * h_i,j * R                         <--- default
# 
#                               + q_d * h_i,j * R ]                       <--- default
# 
# @param short.rate.lattice: interest rate at each node
# @param hazard.rate.lattice: hazard rate at each node
# @param face.value: the face value of the bond
# @param recovery.rate: the recovery rate
# @param n:  number of periods in the binomial model
# @param q:  risk-neutral probability q
# 
# @return bond.price.lattice: n+1 x n+1 matrix of bond prices at each node.
#
buildDefaultableZCBPriceLattice <- function(short.rate.lattice,
                                            hazard.rate.lattice,
                                            face.value,
                                            recovery.rate,
                                            n=ncol(short.rate.lattice)-1,
                                            q=0.5) {


    # n+1 x n+1 matrix (only the upper half is relevant)
    bin.price.lattice <- matrix(NA, nrow=n+1, ncol=n+1)

    # set final prices 
    final.price.vec <- rep(face.value, n+1)
    bin.price.lattice[,n+1] = final.price.vec

    # recovery
    recovery.value <- recovery.rate * face.value

    for (t in n:1) {            # iterating backwards over cols (time periods)
        for (j in 1:t) {        # iterating over rows (states of the world)
            upmove <- bin.price.lattice[j,t+1]
            downmove <- bin.price.lattice[j+1,t+1]

            r <- short.rate.lattice[j,t]
            h <- hazard.rate.lattice[j,t]

            bin.price.lattice[j,t] = 1 / (1 + r) * (q * (1 - h) * upmove + 
                                                    (1 - q) * (1 - h) * downmove + 
                                                    q * h * recovery.value + 
                                                    (1 - q) * h * recovery.value)
        }
    }

    bin.price.lattice

}


#
# 
# A bond price lattice is a binomial lattice where the final prices (t=n) are equal
# to the face value of the bond + the coupon payment.
#
# Coupon payments are NOT discounted, since we assume they are paid in the same
# time period in which they occur.
#
# 
# @param short.rate.lattice: interest rate at each node
# @param face.value: the face value of the bond
# @param coupon.rate: rate of return of the bond (coupons = face.value * coupon.rate )
# @param n:  number of periods in the binomial model
# @param q:  risk-neutral probability q
# 
# @return bond.price.lattice: n+1 x n+1 matrix of bond prices at each node.
#
buildBondPriceLattice <- function(short.rate.lattice,
                                  face.value,
                                  coupon.rate,
                                  n,
                                  q=0.5) {

    buildBinomialPriceLattice( final.price.vec = rep(face.value, n+1 ),
                               short.rate.lattice = short.rate.lattice,
                               coupon.payment.lattice = matrix( face.value * coupon.rate, nrow=n+1, ncol=n+1 ),
                               n = n,
                               q = q)
}



#
# Build a binomial price lattice.
#
# The final prices in the binomial lattice are given by final.price.vec (plus 
# the final coupon payments, if any).
#
# The binomial price lattice is filled in by working backwards thru the lattice, 
# computing prices at each node based on the discounted expected value (under 
# Q probability space) of the asset in the next time period.
#
#
#            [  S_t+1  ]         
# S_t = E_Q_t[ ------- ] 
#            [ 1 + r_t ]
# 
#
# Coupon payments are added to the price at each node.  Coupons are discounted
# if discount.coupons=T (e.g. when coupons are paid in arrears).
#
#
#            [  S_t+1  ]     C_t        <--- coupon for node t
# S_t = E_Q_t[ ------- ] + ---------
#            [ 1 + r_t ]    1 + r_t     <--- discounted if necessary
# 
#
#
# Q probability space is defined by q (and 1-q).  Therefore the price at each
# node is calculated by:
#
#           1                                            C_t
# S_t = ----------- (q * upmove + (1-q) * downmove) + ----------
#        1 + r_t                                       1 + r_t
#
#
# In order to build the price lattice, the user must specify:
#
#   1. the final value of the asset for all nodes at t=n
#       e.g. a bond's final value is its face value (+coupon)
#       e.g. an option's final value is the diff between the underlying and the strike
#   2. the interest rate at every node
#   3. the risk-neutral probabilities associated with upmoves and downmoves in the lattice
#   4. coupon payments (if any) at every node
#   5. whether or not to discount coupon payments (paid in arrears)
#
#
# @param final.price.vec: the final prices (t=n) for the lattice
# @param short.rate.lattice: interest rates at all nodes
# @param n: number of time periods in the lattice
# @param q: risk-neutral probabilities
# @param coupon.payment.lattice: coupon payments for each node
# @param discount.coupons: T or F, whether or not to discount coupon payments
#
# @return bin.price.lattice: n+1 x n+1 matrix
#
#         Note that the j indices (states of the world) are reversed:
#
#           bin.price.lattice[j,t] corresponds to p_t-1,t-j
#
#           bin.price.lattice[1,1] corresponds to p_0,0
#
#           bin.price.lattice[1,2] corresponds to p_1,1
#           bin.price.lattice[2,2] corresponds to p_1,0
#
buildBinomialPriceLattice <- function( final.price.vec = rep(0,n+1),
                                       short.rate.lattice = matrix(0,nrow=n+1,ncol=n+1),
                                       n=length(underlying.final.price)-1, 
                                       q=0.5,
                                       coupon.payment.lattice = matrix(0,nrow=n+1,ncol=n+1), 
                                       discount.coupons=F ) {              

    # discount coupons if necessary
    coupons <- if (discount.coupons) {
        coupon.payment.lattice / (1 + short.rate.lattice[1:(n+1),1:(n+1)]) 
    } else {
        coupon.payment.lattice
    }

    # n+1 x n+1 matrix (only the upper half is relevant)
    bin.price.lattice <- matrix(NA, nrow=n+1, ncol=n+1)

    # set final prices (including coupons)
    bin.price.lattice[,n+1] = final.price.vec[1:(n+1)] + coupons[,n+1]

    for (t in n:1) {            # iterating backwards over cols (time periods)
        for (j in 1:t) {        # iterating over rows (states of the world)
            upmove <- bin.price.lattice[j,t+1]
            downmove <- bin.price.lattice[j+1,t+1]

            r <- short.rate.lattice[j,t]

            bin.price.lattice[j,t] = 1 / (1 + r) * (q * upmove + (1 - q) * downmove) + coupons[j,t]
        }
    }

    bin.price.lattice
}


#
#
# Build a defaultable bond payment schedule.
#
#
#  bond.value = SUM_k=1..n c * q(t_k) * Z[t_k]_0 
# 
#               + F * q(t_n) * Z[t_n]_0
# 
#               + R * F * SUM_k=1..n (q(t_k-1) - q(t)) * Z[t_k]_0
# 
#                     Note: Z[t_k]_0 = d(0,t_k)
#
# @param face.value: face value of bond
# @param coupon.rate: of bond
# @param delta: coupon payment interval (fraction of a year) (default=1/2)
# @param recovery.rate: fraction of face value that is recoverable upon default
# @param hazard.rate.vec: hazard rates for all time periods
# @param n: number of time periods 
# @param r: fixed, deterministic interest rate for all time periods
# @param discount.rate.vec: vector of discount rates for all time periods.
#                           default: 1 / (1 + r*delta)^(0:n),
#
# @return list(hazard.rate = hazard.rate.vec[1:n+1],
#              survival.probability = q,
#              default.probability = default.prob.vec,
#              discount.rate = discount.rate.vec[1:n+1],
#              coupon.payment = coupon.payment.vec,
#              principal.payment = face.value.payment.vec,
#              recovery.payment = recovery.payment.vec,
#              all.payment = coupon.payment.vec + face.value.payment.vec + recovery.payment.vec)
#
buildDefaultableBondSchedule <- function( face.value,
                                          coupon.rate,
                                          delta = 1/2,
                                          recovery.rate,
                                          hazard.rate.vec,
                                          n = length(hazard.rate.vec) - 1,
                                          r,
                                          discount.rate.vec = 1 / (1 + r*delta)^(0:n) ) {

    #
    # survival probabilities:
    # q(t) = q(t-1) * (1 - h_t-1)
    #      = PROD_k=1..n (1 - h_k-1)
    #
    q = c(1, cumprod(1 - hazard.rate.vec[1:n]))

    #
    # default probabilities = diff in survival probabilities
    # dp(t) = q(t-1) - q(t)
    #       = q(t-1) * h_t-1
    #
    # Note: q[1:n] = q[t=0:t=n-1]
    # Note: q[-1] = q[2:end] = q[t=1:t=n]
    # 
    default.prob.vec <- c(0, q[1:n] - q[-1])

    coupon.payment.vec <- rep(face.value * coupon.rate, n+1) * q * discount.rate.vec
    coupon.payment.vec[1] = 0

    face.value.payment <- face.value * q[n+1] * discount.rate.vec[n+1]
    face.value.payment.vec <- c( rep(0,n), face.value.payment)

    recovery.payment.vec <- recovery.rate * face.value * default.prob.vec * discount.rate.vec

    list(hazard.rate = hazard.rate.vec[1:(n+1)],
         survival.probability = q,
         default.probability = default.prob.vec,
         discount.rate = discount.rate.vec[1:(n+1)],
         coupon.payment = coupon.payment.vec,
         principal.payment = face.value.payment.vec,
         recovery.payment = recovery.payment.vec,
         all.payment = coupon.payment.vec + face.value.payment.vec + recovery.payment.vec)

}

#
#
# Compute the price of a defaultable bond.
#
#  bond.value = SUM_k=1..n c * q(t_k) * Z[t_k]_0 
# 
#               + F * q(t_n) * Z[t_n]_0
# 
#               + R * F * SUM_k=1..n (q(t_k-1) - q(t)) * Z[t_k]_0
# 
#                     Note: Z[t_k]_0 = d(0,t_k)
#
# @param face.value: face value of bond
# @param coupon.rate: of bond
# @param delta: coupon payment interval (fraction of a year) (default=1/2)
# @param recovery.rate: fraction of face value that is recoverable upon default
# @param hazard.rate.vec: hazard rates for all time periods
# @param n: number of time periods 
# @param r: fixed, deterministic interest rate for all time periods
# @param discount.rate.vec: vector of discount rates for all time periods.
#                           default: 1 / (1 + r*delta)^(0:n),
#
# @return list(face.value = face.value,
#              coupon.rate = coupon.rate,
#              bond.schedule = bond.schedule,
#              n=n,
#              price=sum(bond.schedule$all.payment))
#
computeDefaultableBondPricing <- function( face.value,
                                           coupon.rate,
                                           delta = 1/2,
                                           recovery.rate,
                                           hazard.rate.vec,
                                           n = length(hazard.rate.vec) - 1,
                                           r,
                                           discount.rate.vec = 1 / (1 + r*delta)^(0:n) ) {

    bond.schedule <- buildDefaultableBondSchedule( face.value = face.value,
                                                   coupon.rate = coupon.rate,
                                                   delta = delta,
                                                   recovery.rate = recovery.rate,
                                                   hazard.rate.vec = hazard.rate.vec,
                                                   n = n,
                                                   r = r,
                                                   discount.rate.vec = discount.rate.vec )

    list(face.value = face.value,
         coupon.rate = coupon.rate,
         bond.schedule = bond.schedule,
         n=n,
         price=sum(bond.schedule$all.payment))
}


#
# Build a CDS pricing schedule
#
# @param N: notional amount
# @param hazard.rate.vec: hazard rates for all time periods
# @param n: number of time periods
# @param delta: premium payment interval (fraction of a year) (default=1/4)
# @param r: fixed, deterministic interest rate for all periods
# @param discount.rate.vec: vector of discount rates for all time periods.
#                           default: 1 / (1 + r*delta)^(0:n),
# @param recovery.rate: percent of notional that is recoverable upon default
#
# @return list(hazard.rate.vec=hazard.rate.vec[1:(n+1)],
#              survival.probability=q,
#              default.probability=default.prob.vec,
#              discount.rate=discount.rate.vec[1:(n+1)],
#              S.par=rep(S.par,n+1), 
#              premium.payment=premiums, 
#              accrued.interest=accrued.interest,
#              protection.payment=protection.payments,
#              premium.leg=premium.leg,
#              protection.leg=protection.leg)
# 
buildCDSSchedule <- function(N,
                             hazard.rate.vec,
                             n = length(hazard.rate.vec) - 1,
                             delta = 1/4,
                             r,
                             discount.rate.vec = 1 / (1 + r*delta)^(0:n),
                             recovery.rate) {

    #
    # survival probabilities:
    # q(t) = q(t-1) * (1 - h_t-1)
    #      = PROD_k=1..n (1 - h_k-1)
    #
    q = c(1, cumprod(1 - hazard.rate.vec[1:n]))

    #
    # default probabilities = diff in survival probabilities
    # dp(t) = q(t-1) - q(t)
    #       = q(t-1) * h_t-1
    #
    # Note: q[1:n] = q[t=0:t=n-1]
    # Note: q[-1] = q[2:end] = q[t=1:t=n]
    # 
    default.prob.vec <- c(0, q[1:n] - q[-1])

    #
    # solving for S_par:
    #
    #            (1 - R) * SUM_k=1..n (q(t_k-1) - q(t_k)) * d(0,t_k)          <--- protection
    # S_par =  ---------------------------------------------------------
    #              delta           
    #              ----- * SUM_k=1..n (q(t_k-1) + q(t_k)) * d(0,t_k)          <--- premiums
    #                2             
    #
    q.sums <- c(0, q[1:n] + q[-1])
    S.par.num <- (1 - recovery.rate) * sum(default.prob.vec * discount.rate.vec) 
    S.par.den <- delta / 2 * sum(q.sums * discount.rate.vec)
    S.par = S.par.num / S.par.den

    #
    # risk-neutral value of premiums: delta * S * N * q(t_k) * d(0,t_k)
    # 
    premiums <- delta * S.par * N * q * discount.rate.vec
    premiums[1] = 0       # no premium paid at t=0

    #
    # accrued interest
    #     delta           
    #  =  ----- * S * N * (q(t_k-1) - q(t_k)) * d(0,t_k)
    #       2             
    # 
    accrued.interest <- delta / 2 * S.par * N * default.prob.vec * discount.rate.vec

    #
    # protection payments
    #   = (1 - R) * N * SUM_k=1..n (q(t_k-1) - q(t_k)) * d(0,t_k)
    #
    protection.payments <- (1 - recovery.rate) * N * default.prob.vec * discount.rate.vec
    
    #
    # If S = S_par, then premium and protection legs should be equal.
    #
    premium.leg <- cumsum(premiums + accrued.interest)
    protection.leg <- cumsum(protection.payments)

    #
    # return everything in a list object
    #
    list(hazard.rate.vec=hazard.rate.vec[1:(n+1)],
         survival.probability=q,
         default.probability=default.prob.vec,
         discount.rate=discount.rate.vec[1:(n+1)],
         S.par=rep(S.par,n+1), 
         premium.payment=premiums, 
         accrued.interest=accrued.interest,
         protection.payment=protection.payments,
         premium.leg=premium.leg,
         protection.leg=protection.leg)
}



#
# Compute monthly mortgage payment
# 
#          c * (1 + c)^n * M_0
#     B = ----------------------
#             (1 + c)^n - 1
#
# @param mortgage.balance (M_0): initial (or remaining) mortgage balance
# @param monthly.mortgage.rate (c):  monthly mortgage coupon rate
# @param n:  number of months remaining in mortgage
#
# @return amortized monthly mortgage payment
#
computeMonthlyMortgagePayment <- function(mortgage.balance,
                                          monthly.mortgage.rate,
                                          n) {
    M0 = mortgage.balance
    c = monthly.mortgage.rate

    B.num <-  c * (1 + c)^n * M0
    B.den <-  (1 + c)^n - 1
    B.num / B.den
}


#
# Compute present value of mortgage
#
#       c * (1 + c)^n * M_0      (1 + r)^n - 1
# PV = ---------------------- * ----------------
#          (1 + c)^n - 1           r * (1 + r)^n
# 
# @return present value of the mortgage, based on market rate
#
computeMortgagePresentValue <- function( monthly.payment,
                                         monthly.market.rate,
                                         n ) {
    B <- monthly.payment
    r <- monthly.market.rate

    B * ((1 + r)^n - 1) / r / (1 + r)^n
}


#
# 
# Assumes: no prepayment, no default
#
#           c * M_0              (1 + r)^n - (1 + c)^n
# PV  = ------------------  *   ------------------------ 
#        (1 + c)^n - 1            (r - c) * (1 + r)^n     
#
#
#
# @return the present value of the mortgage's principal-only payments
#
computeMortgagePrincipalPresentValue <- function( mortgage.balance,
                                                  monthly.mortgage.rate,
                                                  monthly.market.rate,
                                                  n ) {
    c = monthly.mortgage.rate
    r = monthly.market.rate
    M0 = mortgage.balance

    c * M0 * ( (1 + r)^n - (1 + c)^n ) / ( (1 + c)^n - 1 ) / ( (r - c) * (1 + r)^n ) 
}


#
#
# Build a mortgage payment schedule for a level-payment mortgage.
#
# @param mortgage.balance: initial balance
# @param monthly.mortgage.rate: monthly coupon rate for the mortgage
# @param n: number of time periods (months) in the mortgage
#
# @return list(begin.balance=begin.balance.vec,
#              monthly.payment=monthly.payment.vec,
#              mortgage.interest=mortgage.interest.vec,
#              principal.payment=principal.payment.vec,
#              end.balance=end.balance.vec) 
#
buildMortgageSchedule <- function(mortgage.balance,
                                  monthly.mortgage.rate,
                                  n) {
        
    #
    # Monthly payment
    # 
    B <- computeMonthlyMortgagePayment(mortgage.balance = mortgage.balance,
                                       monthly.mortgage.rate = monthly.mortgage.rate,
                                       n = n)
    monthly.payment.vec <- rep(B,n)

    # 
    # Mortgage balance remaining (M_k)
    #
    #                  (1 + c)^n - (1 + c)^k
    #    M_k = M_0 * -------------------------
    #                      (1 + c)^n - 1
    #
    M0 = mortgage.balance
    c = monthly.mortgage.rate
    k <- 1:n

    mortgage.balance.num <- (1 + c)^n - (1 + c)^k
    mortgage.balance.den <- (1 + c)^n - 1

    end.balance.vec <- M0 * mortgage.balance.num / mortgage.balance.den
    begin.balance.vec <- c(M0, end.balance.vec[1:n-1])

    # 
    # Interest and principal
    # 
    # I_k = c * M_k-1         <--- interest payment in k'th period
    #
    # P_k = B - I_k           <--- principal payment in k'th period
    #
    mortgage.interest.vec <- begin.balance.vec * c
    principal.payment.vec <- monthly.payment.vec - mortgage.interest.vec

    return( list(begin.balance=begin.balance.vec,
                 monthly.payment=monthly.payment.vec,
                 mortgage.interest=mortgage.interest.vec,
                 principal.payment=principal.payment.vec,
                 end.balance=end.balance.vec) )
}


#
# Build a Conditional Prepayment Rate schedule using the given parameters.
#
# @param psa.multiplier: multiple of base PSA (default: base PSA = 100)
# @param seasoning: number of time periods that have already passed (default: 0)
# @param n: number of time periods total (including seasoning)
#
# @return cpr.rate[(seasoning+1):n]
#
buildCPRSchedule <- function(psa.multiplier = 100,
                             seasoning = 0,
                             n = 30*12) {

    cpr.rate.vec <- 6 * 1:30/30 
    cpr.rate.vec <- c( cpr.rate.vec, rep(6, max(30*12,n)) )

    cpr.rate.vec <- cpr.rate.vec * psa.multiplier / 100 / 100

    cpr.rate.vec[(seasoning+1):n]
}


#
# Build a Pass-Thru MBS payment schedule
#
# @param mortgage.balance: initial balance of the mortgage pool
# @param monthly.mortgage.rate: monthly mortgage coupon rate
# @param monthly.passthru.rate: monthly passthru rate (paid to MBS investors)
# @param cpr.rate.vec: CPR rate for all time periods
# @param n: number of time periods (months) in the mortgage
#
# @return list(cpr.rate=cpr.rate.vec,
#              smm.rate=smm.rate.vec,
#              begin.balance=begin.balance.vec,
#              monthly.payment=monthly.payment.vec,
#              mortgage.interest=mortgage.interest.vec,
#              passthru.interest=passthru.interest.vec,
#              principal.payment=principal.payment.vec,
#              principal.prepayment=principal.prepayment.vec,
#              total.principal=total.principal.vec,
#              end.balance=end.balance.vec) )
#
buildMBSPassThruSchedule <- function( mortgage.balance,
                                      monthly.mortgage.rate,  
                                      monthly.passthru.rate, 
                                      cpr.rate.vec,   
                                      n = length(cpr.rate.vec) ) {

    smm.rate.vec <- 1 - (1 - cpr.rate.vec)^(1/12)

    begin.balance.vec <- rep(NA,n)
    end.balance.vec <- rep(NA,n)
    monthly.payment.vec <- rep(NA,n)
    mortgage.interest.vec <- rep(NA,n)
    passthru.interest.vec <- rep(NA,n)
    principal.payment.vec <- rep(NA,n)
    principal.prepayment.vec <- rep(NA,n)
    total.principal.vec <- rep(NA,n)

    for (i in 1:n) {

        begin.balance.vec[i] <- mortgage.balance

        monthly.payment.vec[i] <- computeMonthlyMortgagePayment(mortgage.balance = begin.balance.vec[i],
                                                                monthly.mortgage.rate = monthly.mortgage.rate,
                                                                n = n+1-i)

        mortgage.interest.vec[i] <- begin.balance.vec[i] * monthly.mortgage.rate
        passthru.interest.vec[i] <- begin.balance.vec[i] * monthly.passthru.rate

        principal.payment.vec[i] <- monthly.payment.vec[i] - mortgage.interest.vec[i]

        principal.prepayment.vec[i] <- (begin.balance.vec[i] - principal.payment.vec[i]) * smm.rate.vec[i]
        total.principal.vec[i] <- principal.payment.vec[i] + principal.prepayment.vec[i]

        end.balance.vec[i] <- begin.balance.vec[i] - total.principal.vec[i]

        # note: mortgage.balance updated for every iteration
        mortgage.balance <- end.balance.vec[i]
    }

    return( list(cpr.rate=cpr.rate.vec,
                 smm.rate=smm.rate.vec,
                 begin.balance=begin.balance.vec,
                 monthly.payment=monthly.payment.vec,
                 mortgage.interest=mortgage.interest.vec,
                 passthru.interest=passthru.interest.vec,
                 principal.payment=principal.payment.vec,
                 principal.prepayment=principal.prepayment.vec,
                 total.principal=total.principal.vec,
                 end.balance=end.balance.vec) )
}


#
#
# Build a payment schedule for the tranches of a pass-thru MBS
#
# @param mbs: pass-thru MBS payment schedule, as built by buildMBSPassThruSchedule
# @param tranche.begin.belances: initial balances for all tranches.
#                                should sum to the full initial mortgage balance.
# @param n: number of time periods (months) in the MBS passthru schedule.
#
# @return: list( list(begin.balance,
#                     total.principal,
#                     passthru.interest,
#                     end.balance )
#          a list of lists, one per tranche
#
buildMBSTrancheSchedule <- function( mbs,
                                     tranche.begin.balances,
                                     n = length(mbs$begin.balance)) {
    
    # 
    # initialize tranches
    # 
    tranches = vector("list", length(tranche.begin.balances))

    for (j in seq_along(tranche.begin.balances)) {

        tranches[[j]] <- list(begin.balance = rep(0, n),
                              total.principal = rep(0, n),
                              passthru.interest = rep(0, n),
                              end.balance = rep(0, n))
    }

    for (i in 1:n) {        # iterate time periods

        begin.balance <- mbs$begin.balance[i]
        passthru.interest <- mbs$passthru.interest[i]
        total.principal <- mbs$total.principal[i]

        for (j in seq_along(tranche.begin.balances)) {    # iterate thru tranches

            tranche <- tranches[[j]]

            tranche$begin.balance[i] <- ifelse( i > 1, 
                                                tranche$end.balance[i-1],
                                                tranche.begin.balances[j] )
            
            tranche$passthru.interest[i] <- (tranche$begin.balance[i] / begin.balance) * passthru.interest
            tranche$total.principal[i] <- min(total.principal, tranche$begin.balance[i])
            tranche$end.balance[i] <- tranche$begin.balance[i] - tranche$total.principal[i]

            # reduce available total.principal for next tranche
            # by the amount paid into this tranche
            total.principal <- total.principal - tranche$total.principal[i]

            tranches[[j]] <- tranche
        }
    }

    tranches
}


#
# Build a margin account schedule for a futures contract.
#
# @param spot.price.vec: spot prices (=futures prices) of the underlying for all time periods.
#                        quoted in dollars.
# @param num.contracts: quantity of the underlying
# @param initial.margin: initial value of the margin account
# @param maintenance.margin: falling below this value triggers a margin call, which 
#                            replenishes the margin account to its initial value.
# @param n: number of time periods (default: length(spot.price.vec))
#
# @return list(spot.price=spot.price.vec,
#              profit=profit.vec,
#              margin.account=margin.account.vec,
#              margin.call=margin.call.vec)
#
buildMarginAccountSchedule <- function( spot.price.vec, 
                                        num.contracts,
                                        initial.margin,
                                        maintenance.margin,
                                        n = length(spot.price.vec) ) {

    profit.vec <- c(0, diff(spot.price.vec) * num.contracts )

    margin.account.vec <- rep(initial.margin, n) 
    margin.call.vec <- rep(0,n)

    for (i in 2:n) {

        margin.account.vec[i] <- margin.account.vec[i-1] + profit.vec[i]

        if (margin.account.vec[i] < maintenance.margin) {
            margin.call.vec[i] <- initial.margin - margin.account.vec[i]
            margin.account.vec[i] <- initial.margin
        }
    }

    list(spot.price=spot.price.vec,
         profit=profit.vec,
         margin.account=margin.account.vec,
         margin.call=margin.call.vec)
}


