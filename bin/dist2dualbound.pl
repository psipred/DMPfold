#!/usr/bin/perl -w

$pthresh = shift @ARGV;

while (<>)
{
    @fields = split(' ');
    $i = $fields[0];
    $j = $fields[1];

    $nbins = @fields - 2;

    if (abs($i - $j) > 4)
    {
	for ($pmax=$k=0; $k<$nbins; $k++)
	{
	    if ($fields[2 + $k] > $pmax)
	    {
		$pmax = $fields[2 + $k];
		$kmax = $k;
	    }
	}

	$psum = $pmax;

	$k1 = $kmax-1;
	$k2 = $kmax+1;

	$lastk1 = $lastk2 = $kmax;

	while ($psum < $pthresh)
	{
	    if ($k1 >= 0)
	    {
		$p1 = $fields[2 + $k1];
	    }
	    else
	    {
		$p1 = 0;
	    }
	    if ($k2 < $nbins)
	    {
		$p2 = $fields[2 + $k2];
	    }
	    else
	    {
		$p2 = 0;
	    }

	    if ($p1 > $p2)
	    {
		$psum += $p1;
		$lastk1 = $k1;
		$k1--;
	    }
	    else
	    {
		$psum += $p2;
		$lastk2 = $k2;
		$k2++;
	    }
	}

	if ($lastk1 == 0)
	{
	    $dmin = 3.5;
	}
	elsif ($lastk1 < 8)
	{
	    $dmin = 4.0 + 0.5 * $lastk1;
	}
	else
	{
	    $dmin = $lastk1;
	}

	if ($lastk2 == 0)
	{
	    $dmax = 4.5;
	}
	elsif ($lastk2 < 8)
	{
	    $dmax = 0.5 * $lastk2 + 4.5;
	}
	else
	{
	    $dmax = $lastk2 + 1.0;
	}

	if ($dmax < $nbins)
	{
	    print "$i $j $dmin $dmax $psum\n";
	}
    }
}
