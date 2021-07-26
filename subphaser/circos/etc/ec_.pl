#!/usr/bin/perl -w
my ($insvg,$outsvg)=@ARGV;
open F,">$outsvg";
open IN,"<$insvg";
while(<IN>){
	if (/<tspan.*?>(\S+)<\/tspan>/){
		my $name="_$1";
		s/<tspan.*?>\S+<\/tspan>/$name/;
	}
	print F $_;
	
}
close IN;
close F;