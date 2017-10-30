#!/usr/local/bin/php
<?php

/********************************************************************
 *
 * options_regression_testing.php
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

########################################################################################################
# global constants

define('MATCH_EQUALITY', 1);
define('MATCH_EQUALITY_ROUNDED', 2);
define('MATCH_RANGE', 3);
define('MATCH_RANGE_ESTIMATE', 4);

define('ROUNDED_RECISION', 6);


########################################################################################################
# function prototype

function validate_result ($results, $expected_results, & $result_string) {

	$validation_result = true;

	if (!is_array($results)) {
		$results = array($results);
	}

	if (count($results) != count($expected_results)) {
		die('invalid test case!');
	}

	# formatting result string
	$index = 0;
	while (list($key, $val) = each($expected_results)) {

		if ($index != 0) {
			$result_string .= PHP_EOL . '    ';
		}

		$actual_value = $results[$index];
		if ($val[0] == MATCH_EQUALITY_ROUNDED) {
			$actual_value = sprintf('%.' . ROUNDED_RECISION . 'f', $actual_value);
		}

		$result_string .= "'$key' => actual: [$actual_value]; expected: ";
		$current_validation = true;
		
		if ($val[0] == MATCH_EQUALITY) {
			$result_string .= '[' . $val[1] . ']';
			if ($val[1] != $actual_value) {
				$current_validation = false;
			}
		} else if ($val[0] == MATCH_EQUALITY_ROUNDED) {
			$result_string .= '[' . $val[1] . ']';
			if ($val[1] != $actual_value) {
				$current_validation = false;
			}
		} else if ($val[0] == MATCH_RANGE) {
			$result_string .= "(${val[0]}, ${val[1]})";
			if ($val[1] > $actual_value || $actual_value > $val[2]) {
				$current_validation = false;
			}
		} else if ($val[0] == MATCH_RANGE_ESTIMATE) {
			$result_string .= "(${val[0]}, ${val[1]})";
			if ($val[1] > $actual_value || $actual_value > $val[2]) {
				$current_validation = false;
			}
		}

		if ($current_validation) {
			$result_string .= '  => V';
		} else {
			$validation_result = false;
			$result_string .= '  => X';
		}

		++ $index;
	}
	
	return $validation_result;
}

########################################################################################################
# regression test cases

$test_cases = array(
	
	array (
		'case_name' => 'vanila call option pricing (black scholes formula)', 
		'function_name' => 'black_scholes_call', 
		'parameters' => array('S' => 50, 'K' => 55, 'r' => 0.05, 'sigma' => 0.2, 'q' => 0.02, 'T' => 1),
		'results' => array (
						'option_value' => array(MATCH_EQUALITY_ROUNDED, '2.594291')
						)
	),

	array (
		'case_name' => 'vanila put option pricing (black scholes formula)', 
		'function_name' => 'black_scholes_put', 
		'parameters' => array('S' => 50, 'K' => 55, 'r' => 0.05, 'sigma' => 0.2, 'q' => 0.02, 'T' => 1),
		'results' => array (
						'option_value' => array(MATCH_EQUALITY_ROUNDED, '5.901976')
						)
	),

	array (
		'case_name' => 'vanila call option implied volatility (black scholes formula)', 
		'function_name' => 'black_scholes_call_implied_vol', 
		'parameters' => array('call' => 2.59429, 'S' => 50, 'K' => 55, 'r' => 0.05, 'q' => 0.02, 'T' => 1),
		'results' => array (
						'implied_vol' => array(MATCH_EQUALITY_ROUNDED, '0.200000')
						)
	),

	array (
		'case_name' => 'forward start option', 
		'function_name' => 'forward_start_call_option', 
		'parameters' => array('S' => 50,'r' => 0.05, 'sigma' => 0.2, 'q' => 0.02, 'T' => 0.5, 'Tp' => 1),
		'results' => array (
						'option_value' => array(MATCH_EQUALITY_ROUNDED, '3.122433')
						)
	),

	array (
		'case_name' => 'chooser option', 
		'function_name' => 'chooser_option', 
		'parameters' => array('S' => 50, 'Kc' => 55, 'Kp' => 55, 'r' => 0.05, 'sigma' => 0.2, 'q' => 0.02, 
			'T' => 0.5, 'Tc' => 1, 'Tp' => 1),
		'results' => array (
						'option_value' => array(MATCH_EQUALITY_ROUNDED, '7.402774')
						)
	),

);


########################################################################################################
# major chunk

$counter = 0;
$case_sucessful = 0;
$case_warning = 0;
$case_failed = 0;

foreach ($test_cases as $test_case) {
	
	$case_name = $test_case['case_name'];
	$function = $test_case['function_name'];
	$parameters = $test_case['parameters'];
	$expected_results = $test_case['results'];

	# call function through parameter names
	$parameters = array_values($parameters);
	$result = call_user_func_array($function, $parameters);

	$result_string = '';
	$validation = validate_result($result, $expected_results, $result_string);

	if ($validation) {
		$test_status = 'pass';
		++ $case_sucessful;
	} else {
		$test_status = 'fail';
		++ $case_failed;
	}

	++ $counter;

	$parameters = implode(', ', $parameters);

	print <<<RESULT

===============================================================================
Case [$counter] - $case_name
-------------------------------------------------------------------------------
Function:	[$function]
Parameters:	$parameters
Results:	
    $result_string
Status:		[$test_status]
===============================================================================

RESULT;

}

# print testing statistics

	print <<<RESULT

>>>> Regression Testing Statistics <<<<
total: $counter pass: $case_sucessful warning: $case_warning fail: $case_failed


RESULT;

?>
