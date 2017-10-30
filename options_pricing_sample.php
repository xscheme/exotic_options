#!/usr/local/bin/php
<?php

/********************************************************************
 *
 * options_pricing_sample.php
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

include('options.inc.php');

$S = 50;
$K = 55;
$r = 0.05;
$sigma = 0.2;
$q = 0.02;
$T = 1;

echo <<<PARAMETERS

=======================
Sample Input Parameters
=======================

\$S (initial stock price):	$S
\$K (the strike price):		$K
\$r (risk-free rate):		$r
\$sigma (volatility):		$sigma
\$q (dividend yield):		$q
\$T (time to maturity):		$T

PARAMETERS;

$option_value = black_scholes_call ($S, $K, $r, $sigma, $q, $T);
$option_value = sprintf('%.6f', $option_value);

echo <<<PARAMETERS

====================
Sample Function Call
====================

\$option_value = black_scholes_call(\$S, \$K, \$r, \$sigma, \$q, \$T);
\$option_value: $option_value 


PARAMETERS;

?>
