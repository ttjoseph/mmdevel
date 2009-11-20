#!/usr/bin/perl

my @values;

while(<>)
{
	if(/=======/)
	{
		$_ = <>; $_ = <>;
		$_ = <> until /=======/;
	} else
	{
		if(/DV\/DL  =\s+(.+)\s+/)
		{
			push @values, $1;
		}
	}
}

for($i=0; $i<= ($#values - 3); $i++)
{
	print $values[$i] . "\n";
}

