#!/usr/bin/perl
use strict;
use warnings;

my $filename = $ARGV[0];
my @Array;
my @BArray;
my @PArray;
my @VArray;

open(my $File, '<:encoding(UTF-8)', $filename)
    || die 'Could not open file $filename\n';

print("File: $filename\n");

while (my $line = <$File>) {
    chomp $line;

    if($line =~ /^\[(\d+)\]\[(\d+)\]\[(\d+)\] is (-?\d+)/) {
        my $x = $1;
        my $y = $2;
        my $z = $3;
        my $layer = $4;
        if($layer eq "-1"){
            $layer = " ";
        }
        $Array[$x][$y] = $layer;
        #print("[$x][$y]: $layer\n");
    }

    if($line =~ /Boundary\[(\d+)\]\[(\d+)\]\[(\d+)\] is (-?\d+)/) {
        my $x = $1;
        my $y = $2;
        my $z = $3;
        my $layer = $4;
        $BArray[$x][$y] = $layer;
        #print("[$x][$y]: $layer\n");
    }
    if($line =~ /ParticleCell\[(\d+)\]\[(\d+)\]\[(\d+)\] is ((-?\d+)|([A-Z]))/) {
        my $x = $1;
        my $y = $2;
        my $z = $3;
        my $layer = $4;
        if($layer eq "V"){
            $layer = " ";
        }
        $PArray[$x][$y] = $layer;
        #print("[$x][$y]: $layer\n");
    }
    if($line =~ /VoxCell\[(\d+)\]\[(\d+)\]\[(\d+)\] is ((-?\d+)|([A-Z]))/) {
        my $x = $1;
        my $y = $2;
        my $z = $3;
        my $layer = $4;
        $VArray[$x][$y] = $layer;
        #print("[$x][$y]: $layer\n");
    }



}

print "@$_\n" for @Array;
print "\n\n";
print "__________________________________________________________________________________________________\n\n";
#print "@$_\n" for @BArray;
#print "\n\n";
print "@$_\n" for @PArray;
print "\n\n";
print "__________________________________________________________________________________________________\n\n";
print "@$_\n" for @VArray;
