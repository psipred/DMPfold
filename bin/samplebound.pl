#!/usr/bin/perl -w

while (<>)
{
    @fields = split(' ');
    $i = $fields[0];
    $j = $fields[1];

    if (abs($i - $j) > 4)
    {
	$pthresh = rand();
	for ($k=0; $k<20; $k++)
	{
	    $pthresh -= $fields[2 + $k];
	    last if ($pthresh <= 0.0);
	}

	if ($k < 19)
	{
	    $p = $fields[2 + $k];
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
	    print "$i $j $d1 $d2 $p\n"; 
	}
    }
}
