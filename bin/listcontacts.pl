#!/usr/bin/perl

$topomin = shift @ARGV;
$minprob = shift @ARGV;

while (<>)
{
    @fields = split;

    if (abs($fields[1] - $fields[0]) >= $topomin)
    {
	$i = $fields[0];
	$j = $fields[1];
	$score = $fields[4];

	if ($score > $minprob)
	{
	    print $_;
	}
    }
}
