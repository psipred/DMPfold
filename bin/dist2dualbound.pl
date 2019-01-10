#!/usr/bin/perl -w

$samp1 = shift @ARGV;
$samp2 = shift @ARGV;

while (<>)
{
    @fields = split(' ');
    $i = $fields[0];
    $j = $fields[1];

    if (abs($i - $j) > 4)
    {
	for ($p=$k=0; $k<19; $k++)
	{
	    $p += $fields[2 + $k];
	    if ($p > $samp1)
	    {
		if ($k == 0)
		{
		    $d1 = 3.5;
		}
		elsif ($k < 8)
		{
		    $d1 = 4.0 + $k * 0.5;
		}
		else
		{
		    $d1 = $k;
		}
		last;
	    }
	}

	$p -= $fields[2 + $k];

	for (; $k<19; $k++)
	{
	    $p += $fields[2 + $k];
	    if ($p > $samp2)
	    {
		if ($k < 8)
		{
		    $d2 = 4.5 + $k * 0.5;
		}
		else
		{
		    $d2 = $k + 1;
		}
		print "$i $j $d1 $d2 $p\n"; 
		last;
	    }
	}
    }
}
