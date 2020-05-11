#!/usr/bin/perl
# -*-Perl-*-
#Last changed 05/05/2009 repaired the printing of the sequence length
#Last changed Time-stamp: <1999-09-29 13:24:23 ivo>
# produce dot bracket notation of an RNA secondary structure from
# Zuker's .ct file
#btw here is the ref for the server http://rna.tbi.univie.ac.at/cgi-bin/RNAfold.cgi
# Modified by Eluemuno Blyden to notate  b file for in-frame base pairings
# The base is determined to be in frame and matched to another base that 
# is in a codon complementary position (ie modulo 3 values are paired 0:2, 1:1 and 2:0)
# This script produces a modified bracket notation file (bf) in which 
# Frame 1 base pairings are represented as [ and ]
# Frame 2 base pairings are represented as { and }
# Frame 3 base pairings are represented as < and >
# The file also contains code for calculating the statistics of these pairings see ct2bFrameStats.pl  

while (<>) {
    @F = split;
# Look for the = sign in the CT header line and  
# entries for different structures later when $F[4] has changed.
    if ((/ = /) && ($#F!=5)) { 
	if (defined($seq)) {
# Print the fasta header and nucleotide base sequence
	    print ">$seqid bf format $frame0 [1] $frame1 {2} $frame2 <3>\n";
	    print "$seq\n" unless ($seq eq $oldseq);
	    $oldseq = $seq;
# Add up all the base pairs in all frames and print bunch of stats
#SeqID	Energy 	F1count	F2count	F3count	F1/SF	F2/SF	F3/SF	F1/U	F2/U	F3/U	U/Len	U	Len
#	    $sumframes = $frame0 + $frame1 + $frame2;
# Print the bracket notation sequence with the (mfe) at the end
	    print "$s ($E)\n";
#	    print "$seqid\t";
#	    print "$E\t";
#	    print "$frame0\t";
#	    print "$frame1\t";
#	    print "$frame2\t";
# What proportion of the paired bases are in each of the three frames?
#	    $pf0 = $frame0/$sumframes; print "$pf0\t";
#	    $pf1 = $frame1/$sumframes; print "$pf1\t";
#	    $pf2 = $frame2/$sumframes; print "$pf2\t";
# What frame accounts for gains or losses of base pairs during evolution?
#	   $uf0 = $frame0/$unpaired; print "$uf0\t";
#	   $uf1 = $frame1/$unpaired; print "$uf1\t";
#	   $uf2 = $frame2/$unpaired; print "$uf2\t";
#What proportion of bases are unpaired in the whole sequence?
#	    $pu = $unpaired/$F[0]; print "$pu\t";
#	    print "$unpaired\t";
#	    print "$F[0]\t";
# What is the average mfe per bp for the sequence?
#	     print "$N\t";
#	    $Eb = $E/$F[0]; print "$Eb\n";	
	}
#set the seqID, seq length and energy here so they can be printed in table with the right row
	$N = $F[0];	   
	$seqid = $F[4];
	$E = $F[3];
	$s = $seq = "";
	$frame0 = 0;
	$frame1 = 0;
	$frame2 = 0;
	$unpaired = 0;
	next;
    }
    $seq .= $F[1];
    if ($F[4]==0) {$s .= "."; ++$unpaired; next}
    if (($F[0]<$F[4]) and ($F[4] % 3 == 0) and ($F[0] % 3 == 2)) {$s .= "["; ++$frame0; next};
    if (($F[0]<$F[4]) and ($F[4] % 3 == 1) and ($F[0] % 3 == 1)) {$s .= "["; ++$frame0; next};
    if (($F[0]<$F[4]) and ($F[4] % 3 == 2) and ($F[0] % 3 == 0)) {$s .= "["; ++$frame0; next};

    if (($F[0]<$F[4]) and ($F[4] % 3 == 0) and ($F[0] % 3 == 1)) {$s .= "{"; ++$frame1; next};
    if (($F[0]<$F[4]) and ($F[4] % 3 == 1) and ($F[0] % 3 == 0)) {$s .= "{"; ++$frame1; next};
    if (($F[0]<$F[4]) and ($F[4] % 3 == 2) and ($F[0] % 3 == 2)) {$s .= "{"; ++$frame1; next};

    if (($F[0]<$F[4]) and ($F[4] % 3 == 0) and ($F[0] % 3 == 0)) {$s .= "<"; ++$frame2; next};
    if (($F[0]<$F[4]) and ($F[4] % 3 == 1) and ($F[0] % 3 == 2)) {$s .= "<"; ++$frame2; next};
    if (($F[0]<$F[4]) and ($F[4] % 3 == 2) and ($F[0] % 3 == 1)) {$s .= "<"; ++$frame2; next};
    
    $s .= "(" if ($F[0]<$F[4]);

    if (($F[0]>$F[4]) and ($F[4] % 3 == 2) and ($F[0] % 3 == 0)) {$s .= "]"; next};
    if (($F[0]>$F[4]) and ($F[4] % 3 == 1) and ($F[0] % 3 == 1)) {$s .= "]"; next};
    if (($F[0]>$F[4]) and ($F[4] % 3 == 0) and ($F[0] % 3 == 2)) {$s .= "]"; next};
    
    if (($F[0]>$F[4]) and ($F[4] % 3 == 2) and ($F[0] % 3 == 2)) {$s .= "}"; next};
    if (($F[0]>$F[4]) and ($F[4] % 3 == 1) and ($F[0] % 3 == 0)) {$s .= "}"; next};
    if (($F[0]>$F[4]) and ($F[4] % 3 == 0) and ($F[0] % 3 == 1)) {$s .= "}"; next};

    if (($F[0]>$F[4]) and ($F[4] % 3 == 2) and ($F[0] % 3 == 1)) {$s .= ">"; next};
    if (($F[0]>$F[4]) and ($F[4] % 3 == 1) and ($F[0] % 3 == 2)) {$s .= ">"; next};
    if (($F[0]>$F[4]) and ($F[4] % 3 == 0) and ($F[0] % 3 == 0)) {$s .= ">"; next};

    $s .= ")" if ($F[0]>$F[4]);
}
print ">$seqid bf format $frame0 [1] $frame1 {2} $frame2 <3>\n";
print "$seq\n" unless ($seq eq $oldseq);
print "$s ($E)\n" if defined($s);
#$sumframes = $frame0 + $frame1 + $frame2;
#print "$seqid\t";
#print "$E\t";
#print "$frame0\t";
#print "$frame1\t";
#print "$frame2\t";
# What proportion of the paired bases are in each of the three frames?
#	    $pf0 = $frame0/$sumframes; print "$pf0\t";
#	    $pf1 = $frame1/$sumframes; print "$pf1\t";
#	    $pf2 = $frame2/$sumframes; print "$pf2\t";
# What frame accounts for gains or losses of base pairs during evolution?
#	   $uf0 = $frame0/$unpaired; print "$uf0\t";
#	   $uf1 = $frame1/$unpaired; print "$uf1\t";
#	   $uf2 = $frame2/$unpaired; print "$uf2\t";
#What proportion of bases are unpaired in the whole sequence?
#	    $pu = $unpaired/$F[0]; print "$pu\t";
#print "$unpaired\t";

# What is the average mfe per bp for the sequence?
#	    $N = $F[0];	    print "$N\t";
#	    $Eb = $E/$F[0]; print "$Eb\n";	
# End of file
