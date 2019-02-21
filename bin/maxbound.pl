#!/usr/bin/perl -w

while (<>)
{
    @fields = split(' ');
    $i = $fields[0];
    $j = $fields[1];

    if (abs($i - $j) > 4 && $fields[2 + 19] < 0.5)
    {
	$maxp = 0;
	for ($k=0; $k<19; $k++)
	{
	    if ($fields[2 + $k] > $maxp)
	    {
		$maxp = $fields[2 + $k];
		if ($k == 0)
		{
		    $d1 = 3.5;
		    $d2 = 4.5;
		}
		elsif ($k < 8)
		{
		    $d1 = 4.0 + $k * 0.5;
		    $d2 = $d1 + 0.5;
		}
		else
		{
		    $d1 = $k;
		    $d2 = $d1 + 1.0;
		}
	    }
	}

	print "$i $j $d1 $d2 $maxp\n"; 
    }
}
