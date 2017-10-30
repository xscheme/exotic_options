<?php

/********************************************************************
 *
 * options.inc.php
 *
 * Copyright 2017 Shu G. & qa-apps.com
 *
 * Exotic options pricing library
 * https://github.com/xscheme/exotic_options
 *
 * Reference: Kerry Back, A Course in Derivative Securities: 
 *            Introduction to Theory and Computation, Springer, 2010.
 * MathPHP: https://github.com/markrogoyski/math-php
 *********************************************************************/

include('MathPHP/Functions/Support.php');
include('MathPHP/Functions/Special.php');

include('MathPHP/Probability/Distribution/Distribution.php');
include('MathPHP/Probability/Distribution/Continuous/ContinuousDistribution.php');
include('MathPHP/Probability/Distribution/Continuous/Continuous.php');
include('MathPHP/Probability/Distribution/Continuous/Normal.php');

include('MathPHP/Probability/Distribution/Multivariate/Normal.php');

include('MathPHP/NumericalAnalysis/RootFinding/NewtonsMethod.php');
include('MathPHP/NumericalAnalysis/RootFinding/Validation.php');

include('MathPHP/LinearAlgebra/Matrix.php');
include('MathPHP/LinearAlgebra/MatrixFactory.php');
include('MathPHP/LinearAlgebra/SquareMatrix.php');

use MathPHP\Probability\Distribution\Continuous;
use MathPHP\NumericalAnalysis\RootFinding;
use MathPHP\Probability\Multivariate;

function black_scholes_call ($S, $K, $r, $sigma, $q, $T) {
	
	#S = initial stock price 
	#K = the strike price 
	#r = risk-free rate
	#sigma = volatility
	#q = dividend yield
	#T = time to maturity

	$BSCall = 0;

	if ($sigma == 0) {
		$BSCall = max(0,exp(-$q * $T) * $S - exp(-$r * $T)* $K);
	} elseif ($T < 0.000000001) {
		$BSCall = max(0, $S - $K);
	} else {
		$sigsrt = $sigma * sqrt($T);
		$d1=(log($S/$K) + ($r - $q + 0.5 * $sigma * $sigma) * $T) / $sigsrt;
		$d2= $d1 - $sigsrt;

		$normal = new Continuous\Normal(0, 1);
		$N1 = $normal->cdf($d1);
		$N2 = $normal->cdf($d2);

		$BSCall = exp(-$q * $T) * $S * $N1 - exp(-$r * $T) * $K * $N2; 
	}

	return $BSCall;
}

function black_scholes_put ($S, $K, $r, $sigma, $q, $T) {
	
	#S = initial stock price 
	#K = the strike price 
	#r = risk-free rate
	#sigma = volatility
	#q = dividend yield
	#T = time to maturity

	$BSPut = 0;

	if ($sigma == 0) {
		$BSPut = max(0, exp(-$r * $T)* $K - exp(-$q * $T) * $S);
	} elseif ($T < 0.000000001) {
		$BSPut = max(0, $K - $S);
	} else {
		$sigsrt = $sigma * sqrt($T);
		$d1=(log($S/$K) + ($r - $q + 0.5 * $sigma * $sigma) * $T) / $sigsrt;
		$d2= $d1 - $sigsrt;

		$normal = new Continuous\Normal(0, 1);
		$N1 = $normal->cdf(-$d1);
		$N2 = $normal->cdf(-$d2);

		$BSPut = exp(-$r * $T) * $K * $N2 - exp(-$q * $T) * $S * $N1; 
	}

	return $BSPut;
}

function black_scholes_call_implied_vol ($call, $S, $K, $r, $q, $T) {

	#call = call option price
	#S = initial stock price 
	#K = the strike price 
	#r = risk-free rate
	#q = dividend yield
	#T = time to maturity

	if ($call > $S * exp(-$q * $T) || 
		$call < $S * exp(- $q * $T) - $K * exp(- $r * $T)) {
		return false;
	} else {

		$fx = function($sigma) use ($call, $S, $K, $r, $q, $T) {
			return Black_Scholes_Call($S, $K, $r, $sigma, $q, $T);
		};

		// Newton's Method
		$args     = [0.5];  // Parameters to pass to callback function (initial guess, other parameters)
		$target   = $call;   // Value of f(x) we a trying to solve for
		$tol      = 0.00001; // Tolerance; how close to the actual solution we would like
		$position = 0;       // Which element in the $args array will be changed; also serves as initial guess. Defaults to 0.
		
		$sigma = RootFinding\NewtonsMethod::solve($fx, $args, $target, $tol, $position); // Solve for x where f(x) = $target

		return $sigma;
	}

}



function european_call_binomial ($S, $K, $r, $sigma, $q, $T, $N) {
	
	#S = Initial stock price; 
	#K = strike price; 
	#r = risk-free rate; 
	#sigma = volatility; 
	#q = dividend yield; 
	#T = time to maturity; 
	#N = number of time periods;

	$dt = $T / $N;
	$df = exp(-$r * $T);

	#Use the Cox-Ross-Rubinstein parameters;
	$u = exp($sigma * sqrt($dt));
	$d = 1 / $u;
	$pu = (exp(($r - $q) * $dt) - $d) / ($u - $d);
	$pd = 1 - $pu;
	$u2 = $u * $u;

	#Cacluate the stock price at the bottom node (the node corresponding to all
	#down moves; 

	$S = $S * pow($d, $N);
	$prob = pow($pd, $N);

	$call = $prob * max($S - $K, 0);

	#Calculate the other N terms in the sum; 
	for ($i=1; $i<=$N; ++$i) {
		$S = $S * $u2;
		$prob = $prob * $pu * ($N - $i + 1) / ($pd * $i);
		$call = $call + $prob * max($S - $K, 0);
	}

	#Sum up all the nodes at the last period, and discount by e^-rT; 

	return $df * $call;

}

function european_put_binomial ($S, $K, $r, $sigma, $q, $T, $N) {
	
	#S = Initial stock price; 
	#K = strike price; 
	#r = risk-free rate; 
	#sigma = volatility; 
	#q = dividend yield; 
	#T = time to maturity; 
	#N = number of time periods;

	$dt = $T / $N;
	$df = exp(-$r * $T);

	#Use the Cox-Ross-Rubinstein parameters;
	$u = exp($sigma * sqrt($dt));
	$d = 1 / $u;
	$pu = (exp(($r - $q) * $dt) - $d) / ($u - $d);
	$pd = 1 - $pu;
	$u2 = $u * $u;

	#Cacluate the stock price at the bottom node (the node corresponding to all
	#down moves; 

	$S = $S * pow($d, $N);
	$prob = pow($pd, $N);

	$put = $prob * max($K - $S, 0);

	#Calculate the other N terms in the sum; 
	for ($i=1; $i<=$N; ++$i) {
		$S = $S * $u2;
		$prob = $prob * $pu * ($N - $i + 1) / ($pd * $i);
		$put = $put + $prob * max($K - $S, 0);
	}

	#Sum up all the nodes at the last period, and discount by e^-rT; 

	return $df * $put;

}


function american_call_binomial ($S, $K, $r, $sigma, $q, $T, $N) {
	
	#S = Initial stock price; 
	#K = strike price; 
	#r = risk-free rate; 
	#sigma = volatility; 
	#q = dividend yield; 
	#T = time to maturity; 
	#N = number of time periods;

	$dt = $T / $N;
	$df = exp(-$r * $dt);

	#Use the Cox-Ross-Rubinstein parameters;
	$u = exp($sigma * sqrt($dt));
	$d = 1 / $u;
	$pu = (exp(($r - $q) * $dt) - $d) / ($u - $d);
	$pd = 1 - $pu;
	$u2 = $u * $u;

	$calls = array();

	# fix the final pay off
	$ST = $S * pow($d, $N);
	for ($i=1; $i<=$N+1; ++$i) {
		$calls[$i] = max($ST - $K, 0);
		$ST = $ST * $u2;
	}

	# backward induction from time N-1
	for ($i=$N-1; $i>=0; --$i) {
		$St = $S * pow($d, $i);
		for ($j=1; $j<=$i+1; ++$j) {
			$calls[$j] = max($St - $K, ($calls[$j+1] * $pu + $calls[$j] * $pd) * $df);
			$St = $St * $u2;
		}
	}

	return $calls[1];
}


function american_put_binomial ($S, $K, $r, $sigma, $q, $T, $N) {

	#S = Initial stock price; 
	#K = strike price; 
	#r = risk-free rate; 
	#sigma = volatility; 
	#q = dividend yield; 
	#T = time to maturity; 
	#N = number of time periods;

	$dt = $T / $N;
	$df = exp(-$r * $dt);

	#Use the Cox-Ross-Rubinstein parameters;
	$u = exp($sigma * sqrt($dt));
	$d = 1 / $u;
	$pu = (exp(($r - $q) * $dt) - $d) / ($u - $d);
	$pd = 1 - $pu;
	$u2 = $u * $u;

	$puts = array();

	# fix the final pay off
	$ST = $S * pow($d, $N);
	for ($i=1; $i<=$N+1; ++$i) {
		$puts[$i] = max($K - $ST, 0);
		$ST = $ST * $u2;
	}

	# backward induction from time N-1
	for ($i=$N-1; $i>=0; --$i) {
		$St = $S * pow($d, $i);
		for ($j=1; $j<=$i+1; ++$j) {
			$puts[$j] = max($K - $St, ($puts[$j+1] * $pu + $puts[$j] * $pd) * $df);
			$St = $St * $u2;
		}
	}

	return $puts[1];

}



function european_call_mc ($S, $K, $r, $sigma, $q, $T, $M) {
	
	#S = Initial stock price; 
	#K = strike price; 
	#r = risk-free rate; 
	#sigma = volatility; 
	#q = dividend yield; 
	#T = time to maturity; 
	#M = number of simulations;

	# generate M normal distruction numbers
	$random_number = array();
	for ($i=0; $i<$M; ++$i) {
		$random_number[] = get_normal_distribution_num(0, 1);
	}

	$CallV = array();

	for ($i=0; $i<$M; ++$i) {
		$logST = log($S) + ($r - $q - 0.5 * $sigma * $sigma) * $T + $sigma * sqrt($T) * $random_number[$i];
		$ST = exp($logST);
		$CallV[$i] = exp(-$r * $T) * max($ST - $K, 0);
	}

	$call = array_sum($CallV) / $M;

	$squared_sum = 0;
	foreach ($CallV as $val) {
		$squared_sum += $val * $val;
	}

	$stderr = sqrt(($squared_sum - $M * $call * $call) / ($M * ($M-1)));

	$result = array($call, $stderr);

	return $result;

}


function forward_start_call_option ($S, $r, $sigma, $q, $T, $Tp) {
	
	$d1 = ($r - $q + 0.5 * $sigma * $sigma) * sqrt($Tp - $T) / $sigma;
	$d2 = $d1 - $sigma * sqrt($Tp - $T);

	$normal = new Continuous\Normal(0, 1);
	$N1 = $normal->cdf($d1);
	$N2 = $normal->cdf($d2);

	$forward_start_call = exp(-$q * $Tp) * $S * $N1 - exp(-$q * $T - $r * ($Tp - $T)) * $S * $N2;

	return $forward_start_call;

}

function chooser_option ($S, $Kc, $Kp, $r, $sigma, $q, $T, $Tc, $Tp) {

	#S = initial stock price 
	#Kc = the strike price of the call option
	#Kp = the strike price of the put option
	#r = risk-free rate
	#q = dividend yield
	#T = date of the choice mad
	#Tc = time to maturity of the call option
	#Tp = time to maturity of the put option

	$fx = function($s) use ($Kc, $Kp, $r, $sigma, $q, $T, $Tc, $Tp) {
		return Black_Scholes_Call($s, $Kc, $r, $sigma, $q, $Tc - $T) 
			- Black_Scholes_Put($s, $Kp, $r, $sigma, $q, $Tp - $T);
	};

	// Newton's Method
	$args     = [$S];  // Parameters to pass to callback function (initial guess, other parameters)
	$target   = 0;   // Value of f(x) we a trying to solve for
	$tol      = 0.00001; // Tolerance; how close to the actual solution we would like
	$position = 0;       // Which element in the $args array will be changed; also serves as initial guess. Defaults to 0.
	
	$Sstar = RootFinding\NewtonsMethod::solve($fx, $args, $target, $tol, $position); // Solve for x where f(x) = $target
	//echo "sstar $Sstar \n";

	$d1 = (log($S/$Sstar) + ($r - $q + 0.5 * $sigma * $sigma) * $T) / ($sigma * sqrt($T));
	$d2 = $d1 - $sigma * sqrt($T);

	$d1c = (log($S/$Kc) + ($r - $q + 0.5 * $sigma * $sigma) * $Tc) / ($sigma * sqrt($Tc));
	$d2c = $d1c - $sigma * sqrt($Tc);

	$d1p = (log($S/$Kp) + ($r - $q + 0.5 * $sigma * $sigma) * $Tp) / ($sigma * sqrt($Tp));
	$d2p = $d1p - $sigma * sqrt($Tp);

	$chooser_value = exp(-$q * $Tc) * $S * multi_normal_cdf($d1, $d1c, sqrt($T/$Tc))
		- exp(-$r * $Tc) * $Kc * multi_normal_cdf($d2, $d2c, sqrt($T/$Tc))
		+ exp(-$r * $Tp) * $Kp * multi_normal_cdf(-$d2, -$d2p, sqrt($T/$Tp)) 
		- exp(-$q * $Tp) * $S  * multi_normal_cdf(-$d1, -$d1p, sqrt($T/$Tp));

	return $chooser_value;

}


function chooser_option_binomial ($S, $Kc, $Kp, $r, $sigma, $q, $T, $Tc, $Tp, $N) {

	#S = initial stock price 
	#Kc = the strike price of the call option
	#Kp = the strike price of the put option
	#r = risk-free rate
	#q = dividend yield
	#T = date of the choice mad
	#Tc = time to maturity of the call option
	#Tp = time to maturity of the put option
	#N = number of time periods;

	$dt = $T / $N;
	$df = exp(-$r * $T);

	#Use the Cox-Ross-Rubinstein parameters;
	$u = exp($sigma * sqrt($dt));
	$d = 1 / $u;
	$pu = (exp(($r - $q) * $dt) - $d) / ($u - $d);
	$pd = 1 - $pu;
	$u2 = $u * $u;

	#Cacluate the stock price at the bottom node (the node corresponding to all
	#down moves; 

	$S = $S * pow($d, $N);
	$prob = pow($pd, $N);

	$call = $prob * max(black_scholes_call ($S, $Kc, $r, $sigma, $q, $Tc - $T), 
		black_scholes_put ($S, $Kp, $r, $sigma, $q, $Tp - $T));

	#Calculate the other N terms in the sum; 
	for ($i=1; $i<=$N; ++$i) {
		$S = $S * $u2;
		$prob = $prob * $pu * ($N - $i + 1) / ($pd * $i);
		$call = $call + $prob * max(black_scholes_call ($S, $Kc, $r, $sigma, $q, $Tc - $T), 
			black_scholes_put ($S, $Kp, $r, $sigma, $q, $Tp - $T));
	}

	#Sum up all the nodes at the last period, and discount by e^-rT; 

	$chooser_value = $df * $call;

	return $chooser_value;

}

function call_on_call ($S, $K, $r, $sigma, $q, $T, $Kp, $Tp) {

	$fx = function($s) use ($K, $r, $sigma, $q, $T, $Kp, $Tp) {
		return Black_Scholes_Call($s, $Kp, $r, $sigma, $q, $Tp - $T) - $K;
	};

	// Newton's Method
	$args     = [$S];  // Parameters to pass to callback function (initial guess, other parameters)
	$target   = 0;   // Value of f(x) we a trying to solve for
	$tol      = 0.00001; // Tolerance; how close to the actual solution we would like
	$position = 0;       // Which element in the $args array will be changed; also serves as initial guess. Defaults to 0.
	
	$Sstar = RootFinding\NewtonsMethod::solve($fx, $args, $target, $tol, $position); // Solve for x where f(x) = $target
	//echo "sstar $Sstar \n";

	$d1 = (log($S/$Sstar) + ($r - $q + 0.5 * $sigma * $sigma) * $T) / ($sigma * sqrt($T));
	$d2 = $d1 - $sigma * sqrt($T);

	$d1p = (log($S/$Kp) + ($r - $q + 0.5 * $sigma * $sigma) * $Tp) / ($sigma * sqrt($Tp));
	$d2p = $d1p - $sigma * sqrt($Tp);

	$N2 = Normal::CDF($d2, 0, 1);

	$normal = new Continuous\Normal(0, 1);
	$N2 = $normal->cdf($d2);

	$call_on_call_value = -exp(-$r*$T) * $K * $N2 + exp(-$q * $Tp) * $S * multi_normal_cdf($d1, $d1p, sqrt($T / $Tp)) 
		- exp(-$r * $Tp) * $Kp * multi_normal_cdf($d2, $d2p, sqrt($T / $Tp));
	
	return $call_on_call_value;

}


function european_basket_call_mc ($S, $K, $r, $cov, $q, $w, $T, $L, $M) {
	
	#S = L-vector of initial stock prices
	#K = strike price
	#r = risk-free rate
	#Cov = L x L matrix of covariances
	#q = L-vector of dividend yields
	#w = L-vector of basket weights
	#T = time to maturity
	#L = number of assets in the basket
	#M = number of simulations

	$basket_calls = array();

	$logST = array();
	$var = array();
	for ($j=0; $j<$L; ++$j) {
		$var[] = $cov[$j][$j];
	} 

	$multiplier = Cholesky($L, $cov);

	# generate L*M normal distruction numbers
	$random_number = array();
	for ($i=0; $i<$M; ++$i) {
		$random_number[$i] = array();
		for ($j=0; $j<$L; ++$j) {
			$random_number[$i][$j] = get_normal_distribution_num(0, 1);
		}
	}

	for ($i=0; $i<$M; ++$i) {

		$basket_price = 0;

		for ($j=0; $j<$L; ++$j) {

			$z = 0;
			for ($k=0; $k<$L; ++$k) {
				$z += $multiplier[$j][$k] * $random_number[$i][$k];
			}

			$logST = log($S[$j]) + ($r - $q[$j] - 0.5 * $var[$j]) * $T + sqrt($T) * $z;
			$ST = exp($logST);

			$basket_price += $w[$j] * $ST;
		}

		$basket_calls[$i] = exp(-$r * $T) * max($basket_price - $K, 0);

	}

	$call = array_sum($basket_calls) / $M;

	$squared_sum = 0;
	foreach ($basket_calls as $val) {
		$squared_sum += $val * $val;
	}

	$stderr = sqrt(($squared_sum - $M * $call * $call) / ($M * ($M-1)));
	$result = array($call, $stderr);

	return $result;

}


function floating_strike_call_mc ($S, $r, $sigma, $q, $T, $Smin, $N, $M) {

	####!!!!!!! not validated#######

	#S = Initial stock price; 
	#r = risk-free rate; 
	#sigma = volatility; 
	#q = dividend yield; 
	#T = time to maturity; 
	#Smin = minimum during previous life of contract;
	#N = number of periods per simulation; 
	#M = number of simulations; 

	$dt = $T / $N;

	# generate M*N normal distruction numbers
	$random_number = array();
	for ($i=0; $i<$M; ++$i) {
		$random_number[$i] = array();
		for ($j=0; $j<$N; ++$j) {
			$random_number[$i][$j] = get_normal_distribution_num(0, 1);
		}
	}

	$CallV = array();

	$nudt = ($r - $q - 0.5 * $sigma * $sigma) * $dt;

	for ($i=0; $i<$M; ++$i) {

		$logSmin = log($Smin);
		$logS = log($S);

		for ($j=0; $j<$N; ++$j) {
			$logS = $logS + $nudt + $sigma * sqrt($dt) * $random_number[$i][$j];
			$logSmin = min($logSmin, $logS);
		}

		$min = exp($logSmin);
		$ST = exp($logS);

		$CallV[$i] = max(0, $ST - $min) * exp(-$r * $T);

	}

	$call = array_sum($CallV) / $M;

	$squared_sum = 0;
	foreach ($CallV as $val) {
		$squared_sum += $val * $val;
	}

	$stderr = sqrt(($squared_sum - $M * $call * $call) / ($M * ($M-1)));
	$result = array($call, $stderr);

	return $result;

}


function american_spread_put_binomial($S, $K, $r, $sigma, $rho, $q, $T, $N) {

	# S = 2-vector of initial stock prices
	# K = strike price
	# r = risk-free rate
	# sigma = 2-vector of volatilities
	# rho = correlation
	# q = 2-vector of dividend yields
	# T = time to maturity
	# N = number of periods in binomial model

	return 0;
}

###############################################
# internally called utility functions


function get_normal_distribution_num ($mean, $stddev) {
	
	$randmax = 9999; 
	$v = 2; 
	while ($v >= 1 || $v == 0) {

		$u1 = (rand(0, $randmax) / $randmax) * 2.0 - 1.0; 
		$u2 = (rand(0, $randmax) / $randmax) * 2.0 - 1.0; 

		$v = $u1 * $u1 + $u2 * $u2; 
	} 
	
	return $mean + $stddev * ($u1 * sqrt( -2.0 * log($v) / $v)); 

}


function multi_normal_cdf ($d1, $d2, $rho) {

	$script = <<<SCRIPT

	mean <- c(0, 0)
	lower <- c(-Inf, -Inf)
	upper <- c($d1, $d2)

	corr <- diag(2)
	corr[lower.tri(corr)] <- $rho
	corr[upper.tri(corr)] <- $rho

	multi_cdf <- pmvnorm(lower, upper, mean, corr)
	attr(multi_cdf,'error') <- NULL
	attr(multi_cdf,'msg') <- NULL
	print(multi_cdf)

SCRIPT;

	$script_path = '/tmp/R_script_' . getmypid() . '.r';
	file_put_contents($script_path, $script);
	
	$result = trim(shell_exec('Rscript ' . $script_path));
	$result = explode(' ', $result);

	unlink($script_path);

	return $result[1];
}

function Cholesky ($L, $cov) {

	# Inputs are L = number of assets
	# Cov = L x L matrix of covariances

	$a = array();
	for ($i=0; $i<$L; ++$i) {
		$a[$i] = array_fill(0, $L, 0);
	}

	for ($i=0; $i<$L; ++$i) {
		$SumSq = 0;
		for ($h=0; $h<$i; ++$h) {
			$SumSq = $SumSq + $a[$i][$h] * $a[$i][$h];
		}

		$a[$i][$i] = sqrt($cov[$i][$i] - $SumSq);
		for ($j=$i; $j<$L; ++$j) {
			$SumPr = 0;
			for ($h=0; $h<$i; ++$h) {
				$SumPr = $SumPr + $a[$i][$h] * $a[$j][$h];
			}
			$a[$j][$i] = ($cov[$i][$j] - $SumPr) / $a[$i][$i];
		}

	}
		
	return $a;
}

?>