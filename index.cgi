#!/usr/bin/env perl

use G;
use POSIX;
use strict;
use Data::Dumper;
use CGI;
use KyotoCabinet;

use lib qw{ ./lib };
#use TALEN;
use PRIMER;

my $cgi = new CGI;
my $minspacer    = $cgi->param("min")  || 15;
my $maxspacer    = $cgi->param("max")  || 30;
my $minrv        = $cgi->param("rmin") || 15;
my $maxrv        = $cgi->param("rmax") || 20;
my $resultnum    = $cgi->param("resultMax") || 10;
my $jobid        = $cgi->param("jobid");
my $gene         = $cgi->param("gene"); # || "1+zgc:163025";
#my $gene        = $cgi->param("gene") || "6+s1pr5b";
my $sequence     = $cgi->param("sequence");

my $psizemin     = $cgi->param("psizemin") || 20;
my $psizeopt     = $cgi->param("psizeopt") || 22;
my $psizemax     = $cgi->param("psizemax") || 25;
my $primerTM_min = $cgi->param("primerTM_min") || 60;
my $primerTM_opt = $cgi->param("primerTM_opt") || 63;
my $primerTM_max = $cgi->param("primerTM_max") || 65;
my $primerGC_min = $cgi->param("primerGC_min") || 40;
my $primerGC_opt = $cgi->param("primerGC_opt") || 45;
my $primerGC_max = $cgi->param("primerGC_max") || 60;

$resultnum = 100 if ($resultnum > 100);

my $debug = 0;

my $seq;
my ($disnum, $recnum);

#my $fasta = new G("data/s1pr5b_200.fasta", "no cache", "no msg", "fasta");
#$sequence = $fasta->{SEQ};

#$gene = "5+spns3";$sequence = "atgc";

my $weight = 0.9;
my %score; #                                                                    A C G T
$score{"t"} = -log($weight * 57 / 70  + (1.0 - $weight) / 4.0); #NG->T          6 6 1 57  = 70
$score{"c"} = -log($weight * 99 / 107 + (1.0 - $weight) / 4.0); #HD->C          7 99 0 1  = 107
$score{"a"} = -log($weight * 58 / 64  + (1.0 - $weight) / 4.0); #NI->A          58 6 0 0  = 64
$score{"g"} = -log($weight * 26 / 57  + (1.0 - $weight) / 4.0); #NN->G or A     21 8 26 2 = 57

my %mrna;
my $key;
if(!$jobid && !$sequence){
    print $cgi->header();
    &header;
    &form;
}elsif($gene || $sequence || $jobid){
    if($gene && $gene ne "nodata"){
        my ($chr, $name) = ($1, $2) if $gene =~ /(.+)\+(.+)/;
        my $gbk = "./data/dr_ref_Zv9_chr".$chr.".gbk";
        my $gb = load($gbk, "multiple locus", "no msg");
	$seq = substr($gb->{SEQ}, $gb->{$name}->{start} - 101, $gb->{$name}->{end} - $gb->{$name}->{start} + 201);
        $sequence = to_fasta($seq, -name=>$name, -output=>"return");
	$mrna{$gene} = $seq;
    }elsif($sequence){
	$seq = $sequence;
        my $key = $1 if $sequence =~ /(>.+\n)/;
        $seq =~ s/$key//;
        $seq =~ s/[^a-z]+//g;
	$key =~ s/\n//;
	$mrna{$key} = $seq;
    }

    if($jobid){
        print $cgi->header();
        &header;
        &form;
        my $resultfile = "./Results/".$jobid.".xml";
	if(!`grep "length" $resultfile`){
            print qq(
                     <h1>Results</h1>
                     <p>No results.</p>
                     );
        }else{
            print qq(
                     <h1>Results</h1>
                     <div id="grid-example"></div>
                     );
        }
    }else{
        $jobid = sprintf("%010d",time()).sprintf("%05d",$$);
        $!=0;

        my $pid=fork();

        if ( $pid ){
#### Front Job
            print $cgi->header();
            print qq(<meta http-equiv="Refresh" content="1;URL=http://ws.g-language.org/TALEN/running.cgi?jobid=$jobid">);
            close(STDERR);
            close(STDIN);
            close(STDOUT);
            $pid=wait;
        }elsif(defined $pid){



            close(STDERR);
            close(STDIN);
            close(STDOUT);

#### Back Job
	    my $outfile = "./Results/".$jobid.".xml";
	    open(RESULTS, "> $outfile");
	    print RESULTS qq(<?xml version="1.0" encoding="UTF-8"?>
			     <ItemSearchResponse>
			     <Items>
			     );
#####v get over 3000 TALEN results v#####

	    # probability for finding forward TALEN is P(t) * L
            # probability for finding reverse TALEN for each forward TALEN is P(a) * ((maxspacer + maxrv) - (minspacer + minrv))
            # Expect of TALEN pair count is therefore P(t) * L * P(a) * ((maxspacer + maxrv) - (minspacer + minrv))
            # Since limiting factor for TALEN score is P(g), for a given expected count C (default is 3000), 
            # P(g) < C / P(t) * L * P(a) * ((maxspacer + maxrv) - (minspacer + minrv))
	    # we then estimate the necessary cutoff value from this distribution.

	    my ($Pa, $Pt, $Pg, $Pc) = map{$_/length($seq)} seqinfo($seq, -output=>"NULL");
	    my $Ng = ceil($minrv * (3000/($Pa * $Pt * length($seq) * ($maxspacer+$maxrv - $minspacer-$minrv))));
	    my $cutoff = $Ng * $score{"g"} + $Pa * ($minrv - $Ng) * $score{"a"} + $Pt * ($minrv - $Ng) * $score{"t"}  + $Pc * ($minrv - $Ng) * $score{"c"};

	    my %data = &talen_get($seq, $cutoff);

	    warn("RUNNING PRIMER3...\n") if($debug);
#####v get Primer 3 v#####
	    my %primer3_output  = primer3_talen( {'direct'=>$seq, 'comp'=>complement($seq)}, {PRIMER_PRODUCT_SIZE_RANGE=>($minrv * 2 + $minspacer + $psizemin * 2) . '-180', PRIMER_NUM_RETURN=>1000});
	    my %primers_direct  = &parse_primer3($primer3_output{'direct'});
	    my %primers_comp    = &parse_primer3($primer3_output{'comp'}, $seq, 'complement');
	    my %primers = (%primers_direct, %primers_comp);

	    my $tmpfile = "./tmp/tmp_".$jobid.".fastq";
	    my $saifile = "./tmp/tmp_".$jobid.".sai";
	    my $samfile = "./tmp/tmp_".$jobid.".sam";

	    my %fastq;
	    $fastq{"Self"} = substr($seq, 0, 300);
	    for my $id (keys %data){
		$fastq{'TAL' . $data{$id}{TAL1}{dna}} = $data{$id}{TAL1}{dna};
		$fastq{'TAL' . complement($data{$id}{TAL2}{dna})} = complement($data{$id}{TAL2}{dna});
	    }

	    for my $id (keys %primers){
		$fastq{'PRIM' . $primers{$id}{Lseq}} = $primers{$id}{Lseq};
		$fastq{'PRIM' . $primers{$id}{Rseq}} = $primers{$id}{Rseq};
	    }
	    to_fastq(%fastq, -output=>"f", -filename=>$tmpfile);

	    warn("RUNNING BWA...\n") if ($debug);
	    system("/home/t04331nk/bin/bwa-0.6.2/bwa aln ./data/dr_ref.allnuc.fasta $tmpfile -t 8 -n 5 > $saifile ");
	    system("/home/t04331nk/bin/bwa-0.6.2/bwa samse ./data/dr_ref.allnuc.fasta $saifile $tmpfile -n 2000 > $samfile");

	    my $sam = sam_parse($seq, $samfile);

	    say scalar keys %primers if($debug);
	    warn("CHECKING OFF TARGET for PRIMERS...\n") if($debug);
#####v check off target (Primers)  v#####
	    for my $poss (keys %primers){
		my $ouch = 0;
		for my $key ($primers{$poss}{Lseq}, $primers{$poss}{Rseq}){
		    for my $dir (keys %{$sam->{'PRIM' . $key}}){
			for my $chr (keys %{$sam->{'PRIM' . $key}->{$dir}}){
			    for my $query (@{$sam->{'PRIM' . $key}->{$dir}->{$chr}}){
				for my $opponent (@{$sam->{'PRIM' . $primers{$poss}{Lseq}}->{$dir * -1}->{$chr}},
						  @{$sam->{'PRIM' . $primers{$poss}{Rseq}}->{$dir * -1}->{$chr}}){
				    if (abs($query - $opponent) < 10000){
					$ouch=1;
					last;
				    }
				}
				last if($ouch);
			    }
			    last if($ouch);
			}
			last if($ouch);
		    }
		    last if($ouch);
		}
		delete($primers{$poss}) if ($ouch);
	    }
	    say scalar keys %primers if($debug);

	    warn("CHECKING OFF TARGET for TALEN...\n") if($debug);
	    say "talen pairs:",  scalar keys %data if ($debug);
#####v check off target (TALEN) v#####
	    for my $id (keys %data){
		my $ouch = 0;
		for my $key ($data{$id}{TAL1}{dna}, $data{$id}{TAL2}{dna}){
		    for my $dir (keys %{$sam->{'TAL' . $key}}){
			for my $chr (keys %{$sam->{'TAL' . $key}->{$dir}}){
			    for my $query (@{$sam->{'TAL' . $key}->{$dir}->{$chr}}){
				for my $opponent (@{$sam->{'TAL' . $data{$id}{TAL1}{dna}}->{$dir * -1}->{$chr}}, 
						  @{$sam->{'TAL' . $data{$id}{TAL2}{dna}}->{$dir * -1}->{$chr}}){
				    if (abs($query - $opponent) < 500){
					$ouch = 1;
					last;
				    }
				}
				last if($ouch);
			    }
			    last if($ouch);
			}
			last if($ouch);
		    }
		    last if($ouch);
		}
		delete($data{$id}) if($ouch);
	    }
	    say scalar keys %data if ($debug);

	    my %talen_blast;
	    my %bwa_stock;
	    my $num;
	    my $count;


#####v optimization of primer search by limiting the search vector size v#####
	    my @sortedPrimerKeys = sort {$a cmp $b} keys %primers;
	    my @beginPos;
	    my $j = 0;
	    for my $n (0..(int(length($seq)/10))){
		while(1){
		    my ($left, $right, $dir) = split(/_/, $sortedPrimerKeys[$j]);
		    if(floor($left/10) >= $n ){
			last;
		    }else{
			$j ++;
		    }
		    last if ($j >= $#sortedPrimerKeys);
		}
		$beginPos[$n] = $j;
	    }

#####v MAIN PROCESS check for frames and codons with primers v#####
	    for my $id (sort {$data{$a}{TotalScore} <=> $data{$b}{TotalScore}} keys %data){
		$count ++;
		say $count, $num if($count % 1000 == 0);

		my %stock;
		my $blockpos = floor(($data{id}{TAL1}{start} - 100)/100);
		$blockpos = 0 if ($blockpos < 0);
		for my $poss (@sortedPrimerKeys[$beginPos[$blockpos]..$#sortedPrimerKeys]){

		    my $flag_last;
		    my ($left, $right, $dir) = split(/_/, $poss);
		    ($left, $right) = ($primers{$poss}{Lpos}, $primers{$poss}{Rpos});

		    my ($primersLlen, $primersLseq, $primersRlen, $primersRseq) = ($primers{$poss}{Llen}, $primers{$poss}{Lseq}, $primers{$poss}{Rlen}, $primers{$poss}{Rseq});

		    last if ($left > $data{$id}{TAL1}{start});

		    if($left + $primersLlen < $data{$id}{TAL1}{start} && $data{$id}{TAL2}{start} < $right - $primersRlen){
#			next if $dir ne "d"; # YVES

			my $all    = substr($seq, $left, $right - $left + 1);
			my $before = substr($seq, $left, $data{$id}{TAL1}{end} + int(($data{$id}{Spacer}{end} - $data{$id}{Spacer}{start} + 1) / 2) - $left);
			my $after  = substr($seq, $data{$id}{TAL2}{end} - int(($data{$id}{Spacer}{end} - $data{$id}{Spacer}{start} + 1)/2) - 1, $right - ($data{$id}{Spacer}{end} - int(($data{$id}{Spacer}{end} - $data{$id}{Spacer}{start} + 1)/2)) + 1);
			
			my ($llen, $rlen) = ($primersLlen, $primersRlen);
			my ($left2, $right2) = ($left, $right);
			my $method;
			for my $frame (0..5){

			    if($frame % 2 == 1){ # complement
				next unless($dir eq 'c');
			    }else{
				next unless($dir eq 'd');
			    }

			    if($frame != 0 && $frame != 1){
				my $enpri = 0;

				if($frame == 2 || $frame == 3){
				    $enpri = 1;
				}else{
				    $enpri = 2;
				}
				next if $left - $enpri - 1 < 0;
				next if $right + (3 - $enpri) > length($seq);
				$primersLlen = $llen + $enpri;
				$primersLseq = substr($seq, $left - $enpri, $llen + $enpri);
				$primersRlen = $rlen + (3 - $enpri);
				$primersRseq = substr($seq, $right - $rlen + 1, $rlen + (3 - $enpri));
				$all = substr($seq, $left - $enpri, $right - $left + 3 + 1);
# comp $left - $enpri;

				$before = substr($seq, $left - $enpri, $data{$id}{TAL1}{end} + int(($data{$id}{Spacer}{end} - $data{$id}{Spacer}{start} + 1) / 2) - $left + $enpri);
				$after  = substr($seq, $data{$id}{TAL2}{end} - int(($data{$id}{Spacer}{end} - $data{$id}{Spacer}{start} + 1)/2) - 1, $right - ($data{$id}{Spacer}{end} - int(($data{$id}{Spacer}{end} - $data{$id}{Spacer}{start} + 1)/2)) + 1 + (3 - $enpri));
				$left2  = $left - $enpri;
				$right2 = $right + (3 - $enpri);
			    }
			    # direct or complement
			    if($frame % 2 == 1){ # complement
				$all = complement($all);
				my ($after_tmp, $before_tmp) = ($after, $before);
				$after  = complement($before_tmp);
				$before = complement($after_tmp);
			    }

			    # Check no XbaI (TCTAGA) and XpnI (GGTACC) sequences
			    next if $all =~ /(tctaga|ggtacc)/;

			    # Check no at, gt, ta at the end of sequence
			    my $tailseq = substr($all, -2, 2);
			    next if($tailseq  eq "at" || $tailseq eq "gt" || $tailseq  eq "ta");

			    my $afterframe;
			    if(length($after) % 3 == 0){
				$afterframe = 1;
			    }elsif(length($after) % 3 == 1){
				$afterframe = 2;
			    }elsif(length($after) % 3 == 2){
				$afterframe = 0;
			    }

			    if($resultnum =~ /\d+/){
				if(length($all) % 3 == 0){
				    last if $disnum > int($resultnum / 2);
				}else{
				    last if $recnum > int($resultnum / 2);
				}
			    }


			    my $method;
			    while(1){ # check for disruption assay
				last if length($all) % 3 != 0;
				# Check no stop codons after spacer
				last if($after =~ /(taa|tag|tga)/);

				# Check no in-frame start codons after spacer
				last if(substr($after, $afterframe) =~ /(...){0,}(atg|gtg)/);
				
				# Check no in-frame stop codons before spacer
				last if(substr($before, 1) =~ /(...){0,}(taa|tag|tga)/);
				
				$method = "Disruption";
				$disnum ++;
				last;
			    }

			    my $beforeframe;
			    if(length($before) % 3 == 0){
				$beforeframe = 5;
			    }elsif(length($before) % 3 == 1){
				$beforeframe = 3;
			    }elsif(length($before) % 3 == 2){
				$beforeframe = 4;
			    }

			    my $flag2 = 0;
			    while(1){ # check for 1 or 2bp recovery assay
				last if length($all) % 3 == 0;
				# Check no start codons in "after spacer"
				last if $after =~ /(atg|gtg)/;

				# Check no stop codons before spacer, but if there is start codons between the stop-spacer, its OK
				if($before =~ /(tag|tga|taa)(?:...){0,}(atg|gtg){0,}(?:...){0,}.{$beforeframe}$/){
				    last unless(length($2));
				}

				# Check no 1,2-frame (from upstream) stop codons after spacer
				if($after  =~ /(taa|tag|tga)(?:...){0,}(atg|gtg){0,}(?:...){0,}{...}$/){
				    last unless(length($2));
				}

				if(length($all) % 3 == 1){
				    $method = "Recovery (1bp)"; # for 1 bp
				}else{
				    $method = "Recovery (2bp)"; # for 2 bp
				}
				$recnum ++;
				last;
			    }

			    next if !$method;
=head

			    my $rep;
			    $rep .= uc $primersLseq;
			    $rep .= "." x ($data{$id}{TAL1}{start} - ($left2 + $primersLlen) + 1);
			    $rep .= uc $data{$id}{TAL1}{dna};
			    $rep .= $data{$id}{Spacer}{dna};
			    $rep .= uc $data{$id}{TAL2}{dna};
			    $rep .= "." x ($right2 - $data{$id}{TAL2}{start} - $primersRlen - 1);
			    $rep .= uc substr($seq, $right2 - $primersRlen, $primersRlen);
			    
			    next if $stock{$rep};
			    $stock{$rep} ++;
=cut
=head
#####v BWA v#####
			    my $tal = $id."_";
                            my %bwa;
                            my $bwaflag;
                            my $samfile = "./tmp/tmp_".$jobid.".sam";
			    $samfile = './tmp/tmp_134140473724512.sam';
                            for(`grep $tal $samfile`){
                                chomp;
                                my @data = split(/\t/, $_);
                                if(/(TALEN_\d+)_tal1_a/){
                                    my $id = $1;
                                    $bwa{$id}{A}{CHR} = $1 if $data[2] =~ /.+\|(.+)/;
                                    $bwa{$id}{A}{XA} = "Danio_rerio|".$bwa{$id}{A}{CHR}.","."+".$data[3].",1,1\;";
                                    $bwa{$id}{A}{XA} .= $1 if /XA\:Z\:(.+)/;
                                    $bwa{$id}{A}{XA} =~ s/\;/\|tal1\;/g;
                                }elsif(/(TALEN_\d+)_tal2_b/){
                                    my $id = $1;
                                    $bwa{$id}{B}{CHR} = $1 if $data[2] =~ /.+\|(.+)/;
                                    $bwa{$id}{B}{XA} = "Danio_rerio|".$bwa{$id}{B}{CHR}.","."+".$data[3].",1,1\;";
                                    $bwa{$id}{B}{XA} .= $1 if /XA\:Z\:(.+)/;
                                    $bwa{$id}{B}{XA} =~ s/\;/\|tal2\;/g;
                                    
                                    my @tal1 = split(/\;/, $bwa{$id}{A}{XA});
                                    my @tal2 = split(/\;/, $bwa{$id}{B}{XA});
                                    push(@tal1, @tal2);

				    for my $t1 (@tal1){
                                        my ($chr1, $dir1, $pos1, $tal1) = ($1, $2, $3, $4) if $t1 =~ /.+\|(.+).+\,([\-|\+])(\d+).+\|(.+)/;
                                        for my $t2 (@tal1){
                                            my ($chr2, $dir2, $pos2, $tal2) = ($1, $2, $3, $4) if $t2 =~ /.+\|(.+).+\,([\-|\+])(\d+).+\|(.+)/;
                                            if($tal1 eq $tal2){
                                                if($chr1 eq $chr2 && $dir1 ne $dir2){
                                                    if($pos1 - 1000 < $pos2 && $pos1 + 1000 > $pos2){
                                                        $bwaflag ++;
                                                    }
                                                }
                                            }else{
                                                if($chr1 eq $chr2 && $dir1 eq $dir2){
                                                    if($pos1 - 1000 < $pos2 && $pos1 + 1000 > $pos2){
                                                        $bwaflag ++;
                                                    }
                                                }
                                            }
                                            last if $bwaflag > 2;
                                        }
                                        last if $bwaflag > 2;
                                    }
				}
			    }
			    if($bwaflag <= 2){
                                $bwa_stock{$data{$id}{TAL1}{dna}."_".$data{$id}{TAL2}{dna}} = "pass";
                            }else{
                                $bwa_stock{$data{$id}{TAL1}{dna}."_".$data{$id}{TAL2}{dna}} = "hit";
                                last;
                            }

=cut 

			    $num ++;

			    my ($direction, $margeprimer, $tale1_dna_aa, $tale2_dna_aa, $spacer, $tale1_score, $tale2_score);
			    my $all_length = length($all)."bp";
			    my $margeTM = "Forward:_".$primers{$poss}{Ltm}."__Reverse:_".$primers{$poss}{Rtm};
			    my $margeGC = "Forward:_".$primers{$poss}{Lgc}."%__Reverse:_".$primers{$poss}{Rgc}."%";
			    my $margeany = "Forward:_".$primers{$poss}{Lany}."__Reverse:_".$primers{$poss}{Rany};
			    if($frame % 2 == 0){ # direct
				$direction = "direct";
				$margeprimer = "Forward:_".$primersLseq."__Reverse:_".$primersRseq;
#$margeprimer = "Forward:_".$primers{$poss}{Lseq}."__Reverse:_".$primers{$poss}{Rseq};
				$tale1_dna_aa = "DNA:_".$data{$id}{TAL1}{dna}."__Amino Acid:_".seq2rv($data{$id}{TAL1}{dna}, 'direct')."__Position:_".$data{$id}{TAL1}{start}. '..' . $data{$id}{TAL1}{end};
				$tale2_dna_aa = "DNA:_".$data{$id}{TAL2}{dna}."__Amino Acid:_".seq2rv($data{$id}{TAL2}{dna}, 'complement')."__Position:_".$data{$id}{TAL2}{start}. '..' . $data{$id}{TAL2}{end};
				$spacer = 'DNA:_' . $data{$id}{Spacer}{dna} . '__Position:_' . ($data{$id}{TAL1}{end} + 1) . '..' . ($data{$id}{TAL2}{start} - 1);
				$tale1_score = sprintf("%.2f", $data{$id}{TAL1}{score});
				$tale2_score = sprintf("%.2f", $data{$id}{TAL2}{score});
			    }elsif($frame % 2 == 1){
				$direction = "complement";
				$margeprimer = "Forward:_".complement($primersRseq)."__Reverse:_". complement($primersLseq);
#				$margeprimer = "Forward:_".complement($primers{$poss}{Rseq})."__Reverse:_". complement($primers{$poss}{Lseq});
				$tale1_dna_aa = "DNA:_".complement($data{$id}{TAL2}{dna})."__Amino Acid:_".seq2rv($data{$id}{TAL2}{dna}, 'complement')."__Position:_".$data{$id}{TAL2}{start} . '..' . $data{$id}{TAL2}{end};
				$tale2_dna_aa = "DNA:_".complement($data{$id}{TAL1}{dna})."__Amino Acid:_".seq2rv($data{$id}{TAL1}{dna}, 'direct')."__Position:_".$data{$id}{TAL1}{start} . '..' . $data{$id}{TAL1}{end};
				$spacer = 'DNA:_' . complement($data{$id}{Spacer}{dna}) . '__Position:_' . ($data{$id}{TAL2}{end} - 1) . '..' . ($data{$id}{TAL1}{end} + 1);
				$tale1_score = sprintf("%.2f", $data{$id}{TAL2}{score});
				$tale2_score = sprintf("%.2f", $data{$id}{TAL1}{score});
			    }
			    print RESULTS qq(
					     <Item>
					     <ID>set $num</ID>
					     <direction>$direction</direction>
					     <method>$method</method>
					     <Primer>$margeprimer</Primer>
					     <tm>$margeTM</tm>
					     <gc>$margeGC</gc>
					     <any>$margeany</any>
					     <tale1>$tale1_dna_aa</tale1>
					     <tale1score>$tale1_score</tale1score>
					     <spacer>$spacer</spacer>
					     <tale2>$tale2_dna_aa</tale2>
					     <tale2score>$tale2_score</tale2score>
					     <all>$all</all>
					     <length>$all_length</length>
					     </Item>
					     );
			    $flag_last = 1;
			    last;
			}
		    }
		    last if $flag_last;
#		    last if $bwa_stock{$data{$id}{TAL1}{dna}."_".$data{$id}{TAL2}{dna}} eq "hit";
		}
	    if($resultnum =~ /\d+/){
		last if $num == $resultnum;
	    }
 	}
            print RESULTS qq(</Items>
			     </ItemSearchResponse>);
            close RESULTS;
	    
            exit 0;

        }else{
            print "er";
            exit;
        }
        exit;

    }
}else{
    print $cgi->header();
    &header;
    &form;
}

sub header{
    print qq(
             <head>
             <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
             <title>Results</title>
             <link rel="stylesheet" type="text/css" href="./lib/extjs-4.1.0/resources/css/ext-all-gray.css" />
             <link rel="stylesheet" type="text/css" href="./lib/extjs-4.1.0/examples/shared/example.css" />
             <script type="text/javascript" src="./lib/extjs-4.1.0/bootstrap.js"></script>
             
             <script language="javascript">
             var xmlid = $jobid;
             </script>
             
             <!-- page specific -->
             <style type="text/css">
             /* style rows on mouseover */
             .x-grid-row-over .x-grid-cell-inner {
                 font-weight: bold;
             }
             /* shared styles for the ActionColumn icons */
             .x-action-col-cell img {
               height: 16px;
               width: 16px;
               cursor: pointer;
             }
             /* custom icon for the "buy" ActionColumn icon */
             .x-action-col-cell img.buy-col {
                 background-image: url(../shared/icons/fam/accept.gif);
             }
             /* custom icon for the "alert" ActionColumn icon */
             .x-action-col-cell img.alert-col {
                 background-image: url(../shared/icons/fam/error.gif);
             }
             </style>
             <script type="text/javascript" src="./lib/grid-data.js"></script>
             </head>
             );
}


sub talen_get{
    my $upstream  = 't';
    my %data;
    my ($i);
    my $pos = -1;
    my $seq = shift;
    my $cutoff = shift;

#    my %fastq;

    while(1){
	my ($first, $second) = sort {$a<=>$b} (index($seq, 't', $pos + 1), index($seq, 'c', $pos + 1));
	
	if($second < 0){
	    last;
	}elsif($first < 0){
	    $first = $second;
	}
	
	$pos = $first;    
	next if ($pos < $psizemin || $pos > length($seq) - $psizemin - $minrv - $minspacer);
	my $chr = substr($seq, $pos, 1);
	next unless($upstream =~ /$chr/);
	
		
	my $start2 = $pos + $minrv * 2 + $minspacer;
	while(1){
	    my ($first, $second) = sort {$a<=>$b} (index($seq, 'a', $start2 + 1), index($seq, 'g', $start2 + 1));
	    
	    if($second < 0 || $first > $pos + 1 + $maxspacer + 2 * $maxrv){
		last;
	    }elsif($first < 0){
		$first = $second;
	    }
	    
	    last if ($first > length($seq) - $psizemin);

	    $start2 = $first;
	    my $chr2 = substr($seq, $start2, 1);
	    next unless(complement($upstream) =~ /$chr2/);

	    for my $left ($minrv..$maxrv){
		my $talen1 = substr($seq, $pos + 1, $left);
		next if ($talen1 =~ /[^atgc]/);
		my $talen1score = 0;
		$talen1score += $score{$_} for (split(//, lc($talen1)));
                next if ($talen1score > $cutoff);

		for my $right ($minrv..$maxrv){
		    my $center = ($start2 - 1) - ($pos + 1) - $left - $right + 1;
		    next if ($center < $minspacer);
		    my $talen2 = substr($seq, $pos + 1 + $left + $center, $right);
                    my $talen2score = 0;
                    $talen2score += $score{$_} for (split(//, complement(lc($talen2))));
                    next if($talen2score > $cutoff);

		    next if $pos + $left + ceil($center/2) < 30;
		    next if $start2 + 29 > length($seq);
		    $i ++;
		    my $id = "TALEN_".$i;
		    
		    $data{$id}{TAL1}{start} = $pos + 1;
		    $data{$id}{TAL1}{end}   = $pos + $left;
		    $data{$id}{TAL1}{dna}   = $talen1;
		    $data{$id}{TAL1}{score} = $talen1score;
		    
		    $data{$id}{TAL2}{start} = $start2 - 1;
		    $data{$id}{TAL2}{end}   = $pos + 1 + $left + $center;
		    $data{$id}{TAL2}{dna}   = $talen2;
		    $data{$id}{TAL2}{score} = $talen2score;
		    
		    $data{$id}{Spacer}{start} = $pos + $left + 1;
		    $data{$id}{Spacer}{end}   = $pos + $left + $center;
		    $data{$id}{Spacer}{dna}   = substr($seq, $pos, 1) . uc(substr($seq, $pos + 1, $left)) . substr($seq, $pos + 1 + $left, $center) . uc(substr($seq, $pos + 1 + $left + $center, $right)) .  substr($seq, $pos + 1 + $left + $right + $center, 1);
#		    . '_' . uc(substr($seq, $start2-5, 6)) .  '_' . substr($seq, $pos, $left + $center + $right + 3);
		    $data{$id}{length}        = $center;
		    $data{$id}{TotalScore}    = $talen1score + $talen2score;
		    
#		    $fastq{$id."_tal1_a"} = $talen1;
#		    $fastq{$id."_tal2_b"} = complement($talen2);
		}
	    }
	}
    }

    return %data;
}




sub form{
    
    my $seqdata = ">s1pr5b\\n";
    open(SEQDATA, "./data/s1pr5b_200.fasta");
    while(<SEQDATA>){
	chomp;
	next if $_ =~ />/;
	$seqdata .= $_;
    }
    close SEQDATA;

    my $seqdata2 = ">zgc:163025\\n";
    open(SEQDATA2, "./data/zgc163025_200.fasta");
    while(<SEQDATA2>){
	chomp;
	next if $_ =~ />/;
	$seqdata2 .= $_;
    }
    close SEQDATA2;


    my $seqdata3 = ">spns3\\n";
    open(SEQDATA3, "./data/spns3_200.fasta");
    while(<SEQDATA3>){
	chomp;
	next if $_ =~ />/;
	$seqdata3 .= $_;
    }
    close SEQDATA3;


    my $seqdata4 = ">morn3\\n";
    open(SEQDATA4, "./data/morn3_200.fasta");
    while(<SEQDATA4>){
	chomp;
	next if $_ =~ />/;
	$seqdata4 .= $_;
    }
    close SEQDATA4;

    print qq(
             <h1>TALEN, primer getter</h1>
             <form method="POST" action="index.cgi">
             
             <h2>Samples</h2>
	     <input type="button" value="s1pr5b" onclick="sequence.value='$seqdata'; hogefunc();"/>
	     <input type="button" value="zgc:163025" onclick="sequence.value='$seqdata2'; hogefunc();"/>
	     <input type="button" value="spns3" onclick="sequence.value='$seqdata3'; hogefunc();"/>
	     <input type="button" value="morn3" onclick="sequence.value='$seqdata4'; hogefunc();"/>
             <br/>
	     
             );
    
=pod
    <h2>Select gene name</h2>
    my %chr_gene;
    open(GENE, "./data/gene.list");
    while(<GENE>){
        chomp;
        my ($chr, $gene) = split(/\t/, $_);
        my $value = $chr."+".$gene;
	push(@{$chr_gene{$chr}}, $value);
    }
    close GENE;
    for my $chr (sort {$a <=> $b} keys %chr_gene){
	my $i;
	for(@{$chr_gene{$chr}}){
	    $i ++;
	    if($i == 1){
		print qq(
			 <select name="gene">
			 <option value="nodata">== select gene in ==</option>
			 );
	    }elsif($i == scalar(@{$chr_gene{$chr}})){
		print qq(</select>);
	    }else{
		my $name = $_;
		$name =~ s/\+/ /;
		print "<option value=\"".$_."\">".$name."<\/option>\n";
	    }
       }
    }
=cut

    print qq(
             <br/>
             <h2>Submit sequence</h2>
             <textarea name="sequence" id="sequenceid" cols=60 rows=10></textarea>
             <br/>
             <input type="submit" name="submit" value="Submit"/>
             <br/>
             <br/>
             <h2>Options</h2>
             <h3>About TALEN</h3>
             <b>Minimum Spacer Length</b>: <input type="text" size="1" name="min" value="$minspacer"/><br/>
             <b>Maximum Spacer Length</b>: <input type="text" size="1" name="max" value="$maxspacer"/><br/>
             <b>Minimum Repeat Array Length</b>: <input type="text" size="1" name="rmin" value="$minrv"/><br/>
             <b>Maximum Repeat Array Length</b>: <input type="text" size="1" name="rmax" value="$maxrv"/><br/>
             <br/>
             <h3>About Primer</h3>
             <b>Primer Size</b><br />
           Min: <input type="text" name="psizemin" size="1" value="$psizemin"/>
           Opt: <input type="text" name="psizeopt" size="1" value="$psizeopt"/>
           Max: <input type="text" name="psizemax" size="1" value="$psizemax"/>
             <br/><br/>
             <b>Primer Tm</b><br/>
           Min: <input type="text" name="primerTM_min" size="1" value="$primerTM_min"/>
           Opt: <input type="text" name="primerTM_opt" size="1" value="$primerTM_opt"/>
           Max: <input type="text" name="primerTM_max" size="1" value="$primerTM_max"/>
             <br/><br/>
             <b>Primer GC%</b><br/>
           Min: <input type="text" name="primerGC_min" size="1" value="$primerGC_min"/>
           Opt: <input type="text" name="primerGC_opt" size="1" value="$primerGC_opt"/>
           Max: <input type="text" name="primerGC_max" size="1" value="$primerGC_max"/>
             <br/><br/>
	     <b>Number of results</b>: <input type="text" name="resultMax" size="1" value="$resultnum"/>
	     <br/>
<!--	     If you need all results, please input "nolimit". -->
             <br/><br/><br/>
             <input type="submit" name="submit" value="Submit"/>
             );
}

sub seq2score{
    my $seq = shift;
    my $ans;
    for (split(//, lc($seq))){
	$ans += $score{$_};
    }
    return $ans;
}


sub seq2rv{
    my $seq = shift;
    my $direction = shift;
    $seq = complement($seq) if $direction eq "complement";
    $seq =~ s/a/NI /g;
    $seq =~ s/t/NG /g;
    $seq =~ s/c/HD /g;
    $seq =~ s/g/NN /g;

    return $seq;
}

sub parse_primer3{
    my $primer3_output = shift;
    my $seqlen = shift;
    my $direction = shift || 'direct';

    my ($leftpos, $left_seq, $left_tm, $left_gc, $left_any);
    my ($rightpos, $right_seq, $right_tm, $right_gc, $right_any);
    my $no;
    my %primers;
    for my $primerline ( split /\n/, $primer3_output ){
	if($primerline =~ /LEFT PRIMER\s+(\d+).+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+([a-z]+)/){
	    $leftpos = $1;
	    $left_seq = $8; # primer3
	    $left_tm = $3;
	    $left_gc = $4;
	    $left_any = $5;
	}elsif($primerline =~ /RIGHT PRIMER\s+(\d+).+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+([a-z]+)/){
	    $no ++;
	    $rightpos = $1;
	    $right_seq = complement($8); # primer3
	    $right_tm = $3;
	    $right_gc = $4;
	    $right_any = $5;
	
	    if($direction eq 'direct'){
		my $id = sprintf("%07d_%07d_d", $leftpos, $rightpos);
		$primers{$id}{Lseq} = $left_seq;
		$primers{$id}{Rseq} = $right_seq;
		$primers{$id}{Lpos} = $leftpos;
		$primers{$id}{Rpos} = $rightpos;
		$primers{$id}{Llen} = length($left_seq);
		$primers{$id}{Rlen} = length($right_seq);
		
		$primers{$id}{Ltm}  = $left_tm;
		$primers{$id}{Rtm}  = $right_tm;
		$primers{$id}{Lgc}  = $left_gc;
		$primers{$id}{Rgc}  = $right_gc;
		$primers{$id}{Lany} = $left_any;
		$primers{$id}{Rany} = $right_any;
	    }else{
		my $id = sprintf("%07d_%07d_c", $leftpos, $rightpos);
		$primers{$id}{Lseq} = complement($right_seq);
		$primers{$id}{Rseq} = complement($left_seq);
		$primers{$id}{Lpos} = length($seq) - 1 - $rightpos;
		$primers{$id}{Rpos} = length($seq) - 1 - $leftpos;
		$primers{$id}{Llen} = length($right_seq);
		$primers{$id}{Rlen} = length($left_seq);
		
		$primers{$id}{Ltm}  = $right_tm;
		$primers{$id}{Rtm}  = $left_tm;
		$primers{$id}{Lgc}  = $right_gc;
		$primers{$id}{Rgc}  = $left_gc;
		$primers{$id}{Lany} = $right_any;
		$primers{$id}{Rany} = $left_any;
	    }
	}
    }

    return %primers;
}



sub sam_parse{
    my $seq = shift;
    my $sam = shift;
    my %self;
    my $data;
    for my $line (readFile($sam, 1)){
	next if ($line =~ /^\@/);
	
	my @F = split(/\t/, $line);
	
	if($F[0] eq 'Self'){
	    $self{chr} = $F[2];
	    if($F[1] == 16){
		$self{posL} = $F[3] - length($seq);
		$self{posR} = $F[3];
		$self{dir}  = -1;
	    }else{
		$self{posL} = $F[3];
		$self{posR} = $F[3] + length($seq);
		$self{dir}  = 1;
	    }
	}else{
	    my $dir = $F[1] == 16 ? '-' : '+';
	    my $str = $F[2] . ',' . $dir . $F[3] . ',' . $F[5] . '1;';
	    if($line =~ /\tXA:(\S+)/){
		$str .= $1;
	    }
	    
	    for my $seg (split(/;/, $str)){
		my @P = split(/,/, $seg);
		my $pos = substr($P[1], 1);
		my $dir = substr($P[1], 0, 1) eq '+' ? 1 : -1;

		if($P[0] eq $self{chr} && $dir eq $self{dir} && $pos > $self{posL} && $pos < $self{posR}){
#		    say $F[0], $dir, $P[0], $pos, $self{posL}, $self{posR};
		    next;
		}
		
		push(@{$data->{$F[0]}->{$dir}->{$P[0]}}, $pos);
	    }
	}
    }

    return $data;
}

