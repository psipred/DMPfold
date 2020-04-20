#!/usr/bin/perl

# Remove hydrogens and sort/renumber PDB ATOM records according to atom type (DTJ April 2020)

# Extract sort key from ATOM record
sub extractkey
{
    $line = @_[0];

    $resid = substr($line, 22, 5);

    if (substr($line, 13, 3) eq "N  ")
    {
	return "1";
    }
    if (substr($line, 13, 3) eq "CA ")
    {
	return "2";
    }
    if (substr($line, 13, 3) eq "C  ")
    {
	return "3";
    }
    if (substr($line, 13, 3) eq "O  ")
    {
	return "4";
    }
    if (substr($line, 13, 3) eq "OXT")
    {
	return "ZZ";
    }

    $scatm = substr($line, 14, 2);

    $scatm =~ tr/GDEZH/CDEFG/;

    return $scatm;
}

@atomrecs = ();

$lastresid = "";

$atomnum = 1;

while(<>)
{
    if (/^ATOM/)
    {
	if (substr($_, 12, 1) ne "H" && substr($_, 13, 1) ne "H")
	{
	    $resid = substr($_, 22, 5);
	    
	    if ($resid ne $lastresid && @atomrecs > 0)
	    {
		@sorted = sort { extractkey($a) cmp extractkey($b) } @atomrecs;
		foreach $line (@sorted)
		{ 
		    substr($line, 6, 5) = sprintf("%5d", $atomnum);
		    print $line;
		    $atomnum++;
		}
		@atomrecs = ();
	    }
	    
	    $lastresid = $resid;
	    push @atomrecs, $_;
	}
    }
    elsif (/^(TER|END)/)
    {
	@sorted = sort { extractkey($a) cmp extractkey($b) } @atomrecs;
	foreach $line (@sorted)
	{ 
	    substr($line, 6, 5) = sprintf("%5d", $atomnum);
	    print $line;
	    $atomnum++;
	}

	if (/^TER/)
	{
	    print "TER\n";
	}
	else
	{
	    print $_;
	}
	@atomrecs = ();
    }
    else
    {
	print $_;
    }
}
