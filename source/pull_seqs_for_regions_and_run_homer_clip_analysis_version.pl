#!/usr/bin/env perl

use warnings;
use strict;
use Cwd;

# Usage:
#   perl prog BED_FILE [SPECIES] [OUTDIR] [UID}

my $l10p_cutoff = 3;
my $l2fc_cutoff = 3;

my $species = "hg19";
my $output_dir = getcwd();
my $uid = "";
my $window_size = 100;
my $n_bgd_per_event = 5;
# pad 5' end of peaks by n bases
my $fivepad = 0;
my %acceptable_chrs;
for my $n (1..22,"X","Y","M") {
    $acceptable_chrs{"chr".$n} = 1;
}

if (exists $ARGV[1]) {
    $species = $ARGV[1];
    unless ($species eq "hg19") {
	print STDERR "Error - not really tested for non hg19: $species\n";
	exit;
    }
}

if (exists $ARGV[2]) {
    $output_dir = $ARGV[2];
}

$output_dir =~ s/\/$//;
$output_dir .= "/";

unless (-d $output_dir) {
    system("mkdir $output_dir");
}

if (exists $ARGV[3]) {
    $uid = $ARGV[3];
}

my $working_output_dir = $output_dir;
# my $working_output_dir = $output_dir.$uid;
# $working_output_dir =~ s/\/$//;
# $working_output_dir .= "/";
# unless (-d $working_output_dir) {
#     system("mkdir $working_output_dir");
# }

my %all_features;
my %enst2type;
my %enst2ensg;
my %ensg2name;
my %ensg2type;

my $trna_bed = "/storage/vannostrand/genomes/hg19/from_yeolab/hg19-tRNAs.bed";
my $gencode_gtf_file = "/storage/vannostrand/genomes/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf";
my $gencode_tablebrowser_file = "/storage/vannostrand/genomes/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat";
#my $gencode_tablebrowser_file = "/projects/ps-yeolab4/genomes/hg19/gencode_v19/gencodev19_comprehensive";                                                                                               
my $mirbase_fi = "/storage/vannostrand/genomes/hg19/from_yeolab/mirbase.v20.hg19.gff3";
my $lncrna_tablefile = "/storage/vannostrand/genomes/hg19/from_yeolab/lncipedia_5_0_hg19.bed.parsed_ucsc_tableformat";
my $lncrna_fullfi = "/storage/vannostrand/genomes/hg19/from_yeolab/lncipedia_5_0_hg19.gff.parsed";
my $gtf_proteincoding_flag = "no";
if ($species eq "hg19") {
} elsif ($species eq "mm9") {
    $trna_bed = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/mm9-tRNAs.bed";
    $gencode_gtf_file = "/projects/ps-yeolab4/genomes/mm9/gencode.vM1.annotation.gtf";
    $gencode_tablebrowser_file = "/projects/ps-yeolab4/genomes/mm9/gencode.vM1.annotation.gtf.parsed_ucsc_tableformat";
} elsif ($species eq "mm10") {
    $trna_bed = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/mm10-tRNAs.bed";
    $gencode_gtf_file = "/projects/ps-yeolab4/genomes/mm10/gencode/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf";
    $gencode_tablebrowser_file = "/projects/ps-yeolab4/genomes/mm10/gencode/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat";
    $mirbase_fi = "/home/elvannostrand/data/miRBase/mmu_mm10.gff3";
    $lncrna_tablefile = "";
    $lncrna_fullfi = "";
    $gtf_proteincoding_flag = "no";

} elsif ($species eq "hg38") {
    $trna_bed = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/hg38-tRNAs.bed";
    $gencode_gtf_file = "/projects/ps-yeolab4/genomes/GRCh38/gencode/v26/gencode.v26.chr_patch_hapl_scaff.annotation.gtf";
    $gencode_tablebrowser_file = "/projects/ps-yeolab4/genomes/GRCh38/gencode/v26/gencode.v26.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat";
    $mirbase_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/mirbase.v21.hg38.gff3";
    $lncrna_tablefile = "/home/elvannostrand/data/clip/CLIPseq_analysis/lncRNAs/lncipedia_5_0_hg38.bed.parsed_ucsc_tableformat";
    $lncrna_fullfi = "/home/elvannostrand/data/clip/CLIPseq_analysis/lncRNAs/lncipedia_5_0_hg38.gff.parsed";

} elsif ($species eq "rn6") {
    $gencode_gtf_file = "/home/elvannostrand/RN6.gtf";
    $gencode_tablebrowser_file = "/home/elvannostrand/RN6.refGene.table";
    $mirbase_fi = "/home/elvannostrand/rn6_mirbase.gff3";
    $gtf_proteincoding_flag = "all_protein_coding";
    $lncrna_tablefile = "";
    $lncrna_fullfi = "";
    $trna_bed = "";
} else {
    die "species $species not implemented\n";
}
&read_lncrna_parsed($lncrna_fullfi) if ($lncrna_fullfi);
&read_mirbase($mirbase_fi) if ($mirbase_fi);
&parse_trna_list($trna_bed);
&read_gencode_gtf($gencode_gtf_file,$gtf_proteincoding_flag);
&read_gencode($gencode_tablebrowser_file);
&read_gencode($lncrna_tablefile) if ($lncrna_tablefile);

my %annotation_hash;
my %final_feature_list;
my %final_feature_list_bgd;
my $bed_fi = $ARGV[0];
#bed fi is *.annotated_proxdist_miRlncRNA
open(BED,$bed_fi) || die "BED file $bed_fi may not be a file or does not exist.\n";
for my $line (<BED>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    my $chr = $tmp[0];
    my $str = $tmp[5];
    my $start = $tmp[1];
    my $stop = $tmp[2];

    my $l10pval = $tmp[3];
    my $l2fc = $tmp[4];
    if ($l10pval >= $l10p_cutoff && $l2fc >= $l2fc_cutoff) {
    } else {
	next;
    }

## extend peaks 10 bp in the 5' direction
    if ($str eq "+") {
	$start -= $fivepad;
    } elsif ($str eq "-") {
	$stop += $fivepad;
    }
    
    my $final_annot = $tmp[7];
    my ($annotation,$ensg) = split(/\|\|/,$final_annot);

    my $peak_size = $stop - $start;
    
    if ($peak_size > 1000) {
	print STDERR "Warning - peak is greater than 1kb, shouldn't happen with CLIPper peaks - abort\n$bed_fi $line\n";
	next;
    }
    next if ($annotation eq "intergenic");
    $annotation_hash{$annotation} = 1;
    push @{$final_feature_list{$chr}{events}},$start."-".$stop."|".$str;
    push @{$final_feature_list{$chr}{label}},"expt";
    push @{$final_feature_list{$chr}{annotation}},$annotation;
    
#    print "$line\n";
    my $i=0;
    while ($i < $n_bgd_per_event) {
#    for (my $i=0;$i<$n_bgd_per_event;$i++) {
	my $randi = int(rand(scalar(@{$all_features{$annotation}})));
	my $rand_event = $all_features{$annotation}[$randi];
	my ($featureid,$featuretype,$featurechr,$featurepos,$featurestr) = split(/\|/,$rand_event);
	next unless (exists $acceptable_chrs{$featurechr});

	my ($featurestart,$featurestop) = split(/\-/,$featurepos);
	my $featurelen = $featurestop - $featurestart;
	my $rand_position = int(rand($featurelen));

	my $start_offset = int($peak_size / 2);
	my $stop_offset = $peak_size - $start_offset;

	my $final_start = $featurestart + $rand_position - $start_offset;
	my $final_stop = $featurestart + $rand_position + $stop_offset;
	
#	print "bgd $i\t$rand_event\t$final_start\t$final_stop\n";
	push @{$final_feature_list{$featurechr}{events}},$final_start."-".$final_stop."|".$featurestr;
	push @{$final_feature_list{$featurechr}{label}},"bgd";
	push @{$final_feature_list{$featurechr}{annotation}},$annotation;
	$i++;
    }
#    print "\n";

}
close(BED);

my %filehandles;
for my $label ("bgd", "expt") {
    for my $annot ("all",keys %annotation_hash) {
	
	my $outfi = $working_output_dir.$uid.".".$annot.".".$label.".fa";
	open(my $fh, '>', $outfi);
	$filehandles{$label."|".$annot} = $fh;
    }
}

for my $chr (keys %final_feature_list) {
#    print STDERR "doing $chr\n";
    my $chr_fi = "/storage/vannostrand/genomes/".$species."/chromosomes/".$chr.".fa";
    my $chr_seq;
    open(C,$chr_fi) || die "no $chr_fi\n";
    while (<C>) {
	chomp($_);
	next if ($_ =~ /^\>/);
	$chr_seq .= uc($_);
    }
    close(C);

    for (my $j=0;$j<@{$final_feature_list{$chr}{events}};$j++) {
	my $label = $final_feature_list{$chr}{label}[$j];
	my $event = $final_feature_list{$chr}{events}[$j];
	my ($eventpos,$eventstr) = split(/\|/,$event);
	my ($eventstart,$eventstop) = split(/\-/,$eventpos);

	my $annot = $final_feature_list{$chr}{annotation}[$j];

	my $seq = substr($chr_seq,$eventstart,$eventstop - $eventstart);
	if ($eventstr eq "-") {
	    $seq = &revcomp($seq);
	}

	print { $filehandles{$label."|".$annot} } ">".$event."\n".$seq."\n";
	print { $filehandles{$label."|all"} } ">".$event."\n".$seq."\n";
    }
}


for my $label ("bgd", "expt") {
    for my $annot ("all", keys %annotation_hash) {
        close($filehandles{$label."|".$annot});
    }
}
	
for my $annot ("all", keys %annotation_hash) {
    my $outexptfi = $working_output_dir.$uid.".".$annot.".expt.fa";
    my $outbgdfi = $working_output_dir.$uid.".".$annot.".bgd.fa";
    
    my $homer_out_dir = $working_output_dir.$annot."/";
    print STDERR "findMotifs.pl $outexptfi fasta $homer_out_dir -nofacts -p 4 -rna -S 20 -len 5,6,7,8,9 -noconvert -nogo -fastaBg $outbgdfi\n";
#    system("perl /projects/ps-yeolab4/software/eclipconda/envs/eclipanalysis-0.0.3a/bin/findMotifs.pl $outexptfi fasta $homer_out_dir -nofacts -p 4 -rna -S 20 -len 5,6,7,8,9 -noconvert -nogo -fastaBg $outbgdfi");
    system("findMotifs.pl $outexptfi fasta $homer_out_dir -nofacts -p 4 -rna -S 20 -len 5,6,7,8,9 -noconvert -nogo -fastaBg $outbgdfi");

    my $homer_out = $homer_out_dir."homerMotifs.all.motifs";
    system("compareMotifs.pl $homer_out $homer_out_dir");

    my $homer_html = $homer_out_dir."homerResults.html";
    my $new_homer_html = $working_output_dir.$uid.".".$annot.".homerResults.html";
    system("cp $homer_html $new_homer_html");
}

sub revcomp {
    my $seq = shift;
    my $rev = reverse(uc($seq));
    $rev =~ tr/ACGT/TGCA/;
    return($rev);
}

sub parse_trna_list {
    my $trna_fi = shift;
    open(TRNA,$trna_fi) || die "no $trna_fi\n";
    for my $line (<TRNA>) {
        chomp($line);
        my @tmp = split(/\t/,$line);
        my $chr = $tmp[0];
        # trna file is 1 based closed ended - shifting to bed format [0 base open ended) here                                                                                                           
        my $start = $tmp[1]-1;
        my $stop = $tmp[2];
        my $str = $tmp[5];
        my $id = $tmp[3];


        my $feature = $id."|tRNA|".$chr."|".$start."-".$stop."|".$str;
        $enst2ensg{$id}= $id;
        $ensg2name{$id}{$id} = 1;

	push @{$all_features{"tRNA"}},$feature;

    }
    close(TRNA);

}


sub read_mirbase {
    my $mirbase_file = shift;
    open(MIR,$mirbase_fi) || die "no $mirbase_fi\n";
    for my $line (<MIR>) {
        chomp($line);
        $line =~ s/\r//g;
        next if ($line =~ /^\#/);
        my @tmp = split(/\t/,$line);
        if ($tmp[2] eq "miRNA_primary_transcript") {
            my $chr = $tmp[0];
            my $start = $tmp[3]-1;
            my $stop = $tmp[4];
            my $str = $tmp[6];

            if ($tmp[8]=~ /ID\=(\S+?)\;.+Name\=(\S+?)$/ || $tmp[8] =~ /ID\=(\S+?)\;.+Name\=(\S+?)\;/) {
                my $id = $1;
                my $gname = $2;

                my $feature = $id."|miRNA|".$chr."|".$start."-".$stop."|".$str;
#               print STDERR "feature $feature $chr $gname $str\n" if ($gname =~ /mir-21$/);                                                                                                            
		push @{$all_features{"miRNA"}},$feature;

                my $prox_upregion = ($start-500)."|".$start;
                my $prox_dnregion = ($stop)."|".($stop+500);
                for my $proxregion ($prox_upregion,$prox_dnregion) {
                    my ($prox_start,$prox_stop) = split(/\|/,$proxregion);

                    my $prox_feature = $id."|miRNA_proximal|".$chr."|".$prox_start."-".$prox_stop."|".$str;
		    push @{$all_features{"miRNA_proximal"}},$prox_feature;

                }

            } else {
                print STDERR "didn't parse this properly $tmp[8] $line\n";
            }
        }
    }
    close(MIR);
}



sub read_lncrna_parsed {
    my $lncfi = shift;
    open(LN,$lncfi) || die "no $lncfi\n";
    for my $line (<LN>) {
        chomp($line);
        my @tmp = split(/\t/,$line);
        my $enst_id = $tmp[2];
        my $ensg_id = "lncRNA|".$tmp[1];
        
        my $gene_name = $tmp[1];
        my $transcript_type = "lncRNA";
        my $gene_type = "lncRNA";
        
        $enst2ensg{$enst_id} = $ensg_id;
        $ensg2name{$ensg_id}{$gene_name}=1;
        $ensg2type{$ensg_id}{$gene_type}=1;
        $enst2type{$enst_id} = $transcript_type;
        
    }
    close(LN);

}


sub read_gencode {
    ## eric note: this has been tested for off-by-1 issues with ucsc brower table output!                                                                                                             \
                                                                                                                                                                                                       
    my $fi = shift;
#    my $fi = "/projects/ps-yeolab4/genomes/hg19/gencode_v19/gencodev19_comprehensive";                                                                                                                 
    print STDERR "Reading in $fi\n";
    open(F,$fi);
    while (<F>) {
        chomp($_);
        my @tmp = split(/\t/,$_);
        my $enst = $tmp[1];
        next if ($enst eq "name");
        my $chr = $tmp[2];
        my $str = $tmp[3];
        my $txstart = $tmp[4];
        my $txstop = $tmp[5];
        my $cdsstart = $tmp[6];
        my $cdsstop = $tmp[7];

        my @starts = split(/\,/,$tmp[9]);
        my @stops = split(/\,/,$tmp[10]);

        my @tmp_features;

        my $transcript_type = $enst2type{$enst};
        unless ($transcript_type) {
            print STDERR "error transcript_type $transcript_type $enst\n";
        }
        if ($transcript_type eq "protein_coding") {

            for (my $i=0;$i<@starts;$i++) {
                if ($str eq "+") {
                    if ($stops[$i] < $cdsstart) {
                        # exon is all 5' utr
                        push @tmp_features,$enst."|5utr|".$chr."|".($starts[$i])."-".$stops[$i]."|".$str;
                    } elsif ($starts[$i] > $cdsstop) {
                        #exon is all 3' utr                                                                                                                                                                                                       
                        push @tmp_features,$enst."|3utr|".$chr."|".($starts[$i])."-".$stops[$i]."|".$str;
                    } elsif ($starts[$i] > $cdsstart && $stops[$i] < $cdsstop) {
                        #exon is all coding                                                                                                                                                                                                       
                        push @tmp_features,$enst."|CDS|".$chr."|".($starts[$i])."-".$stops[$i]."|".$str;
                    } else {
                        my $cdsregion_start = $starts[$i];
                        my $cdsregion_stop = $stops[$i];

                        if ($starts[$i] <= $cdsstart && $cdsstart <= $stops[$i]) {
                            #cdsstart is in exon                                                                                                                                                                                                       
                            my $five_region = ($starts[$i])."-".$cdsstart;
                            push @tmp_features,$enst."|5utr|".$chr."|".$five_region."|".$str;
                            $cdsregion_start = $cdsstart;
                        }

                        if ($starts[$i] <= $cdsstop && $cdsstop <= $stops[$i]) {
                            #cdsstop is in exon                                                                                                                                                                                                       
                            my $three_region = ($cdsstop)."-".$stops[$i];
                            push @tmp_features,$enst."|3utr|".$chr."|".$three_region."|".$str;
                            $cdsregion_stop = $cdsstop;
                        }

                        my $cds_region = ($cdsregion_start)."-".$cdsregion_stop;
                        push @tmp_features,$enst."|CDS|".$chr."|".$cds_region."|".$str;
                    }
                } elsif ($str eq "-") {
		    if ($stops[$i] < $cdsstart) {
                        # exon is all 5' utr                                                                                                                                                                                                       
                        push @tmp_features,$enst."|3utr|".$chr."|".($starts[$i])."-".$stops[$i]."|".$str;
                    } elsif ($starts[$i] > $cdsstop) {
                        #exon is all 3' utr                                                                                                                                                                                                        
                        push @tmp_features,$enst."|5utr|".$chr."|".($starts[$i])."-".$stops[$i]."|".$str;
                    } elsif ($starts[$i] > $cdsstart && $stops[$i] < $cdsstop) {
                        #exon is all coding                                                                                                                                                                                                        
                        push @tmp_features,$enst."|CDS|".$chr."|".($starts[$i])."-".$stops[$i]."|".$str;
                    } else {
                        my $cdsregion_start = $starts[$i];
                        my $cdsregion_stop = $stops[$i];

                        if ($starts[$i] <= $cdsstart && $cdsstart <= $stops[$i]) {
                            #cdsstart is in exon                                                                                                                                                                                                       
                            my $three_region = ($starts[$i])."-".$cdsstart;
                            push @tmp_features,$enst."|3utr|".$chr."|".$three_region."|".$str."|".$str;
                            $cdsregion_start = $cdsstart;
                        }

                        if ($starts[$i] <= $cdsstop && $cdsstop <= $stops[$i]) {
                            #cdsstop is in exon                                                                                                                                                                                                        
                            my $five_region = ($cdsstop)."-".$stops[$i];
                            push @tmp_features,$enst."|5utr|".$chr."|".$five_region."|".$str."|".$str;
                            $cdsregion_stop = $cdsstop;
                        }

                        my $cds_region = ($cdsregion_start)."-".$cdsregion_stop;
                        push @tmp_features,$enst."|CDS|".$chr."|".$cds_region."|".$str."|".$str;
                    }
                }
            }
            for (my $i=0;$i<scalar(@starts)-1;$i++) {
                # full intron is ($stops[$i]+1)."-".$starts[$i+1]                                                                                                                                      
                # prox is 500bp                                                                                                                                                                        
		
		if ($starts[$i+1]-$stops[$i] > 2 * 500) {
		    if ($str eq "+") {
			push @tmp_features,$enst."|5ss|".$chr."|".($stops[$i])."-".($stops[$i]+$window_size)."|".$str;
			push @tmp_features,$enst."|3ss|".$chr."|".($starts[$i+1]-$window_size)."-".$starts[$i+1]."|".$str;
		    } elsif ($str eq "-") {
			push @tmp_features,$enst."|3ss|".$chr."|".($stops[$i])."-".($stops[$i]+$window_size)."|".$str;
			push @tmp_features,$enst."|5ss|".$chr."|".($starts[$i+1]-$window_size)."-".$starts[$i+1]."|".$str;
		    }
		    
		    
		    push @tmp_features,$enst."|proxintron|".$chr."|".($stops[$i]+$window_size)."-".($stops[$i]+500)."|".$str;
		    push @tmp_features,$enst."|proxintron|".$chr."|".($starts[$i+1]-500)."-".($starts[$i+1]-$window_size)."|".$str;
		    push @tmp_features,$enst."|distintron|".$chr."|".($stops[$i]+500)."-".($starts[$i+1]-500)."|".$str;
		} else {
		    my $midpoint = int(($starts[$i+1]+$stops[$i])/2);

		    if ($starts[$i+1]-$stops[$i] > 2 * $window_size) {
			push @tmp_features,$enst."|proxintron|".$chr."|".($stops[$i]+$window_size)."-".($starts[$i+1]-$window_size)."|".$str;
			
			if ($str eq "+") {
			    push @tmp_features,$enst."|5ss|".$chr."|".($stops[$i])."-".($stops[$i]+$window_size)."|".$str;
			    push @tmp_features,$enst."|3ss|".$chr."|".($starts[$i+1]-$window_size)."-".$starts[$i+1]."|".$str;
			} elsif ($str eq "-") {
			    push @tmp_features,$enst."|3ss|".$chr."|".($stops[$i])."-".($stops[$i]+$window_size)."|".$str;
			    push @tmp_features,$enst."|5ss|".$chr."|".($starts[$i+1]-$window_size)."-".$starts[$i+1]."|".$str;
			}
		    } else {
			if ($str eq "+") {
			    push @tmp_features,$enst."|5ss|".$chr."|".($stops[$i])."-".($midpoint)."|".$str;
			    push @tmp_features,$enst."|3ss|".$chr."|".($midpoint)."-".$starts[$i+1]."|".$str;
			} elsif ($str eq "-") {
			    push @tmp_features,$enst."|3ss|".$chr."|".($stops[$i])."-".($midpoint)."|".$str;
			    push @tmp_features,$enst."|5ss|".$chr."|".($midpoint)."-".$starts[$i+1]."|".$str;
			}
		    }
		}
	    }
	} else {
	    

	    for (my $i=0;$i<@starts;$i++) {
                push @tmp_features,$enst."|noncoding_exon|".$chr."|".($starts[$i])."-".$stops[$i]."|".$str;
            }
            for (my $i=0;$i<scalar(@starts)-1;$i++) {
                if ($starts[$i+1]-$stops[$i] > 2 * 500) {
                    if ($str eq "+") {
                        push @tmp_features,$enst."|noncoding_5ss|".$chr."|".($stops[$i])."-".($stops[$i]+$window_size)."|".$str;
                        push @tmp_features,$enst."|noncoding_3ss|".$chr."|".($starts[$i+1]-$window_size)."-".$starts[$i+1]."|".$str;
                    } elsif ($str eq "-") {
                        push @tmp_features,$enst."|noncoding_3ss|".$chr."|".($stops[$i])."-".($stops[$i]+$window_size)."|".$str;
                        push @tmp_features,$enst."|noncoding_5ss|".$chr."|".($starts[$i+1]-$window_size)."-".$starts[$i+1]."|".$str;
                    }


                    push @tmp_features,$enst."|noncoding_proxintron|".$chr."|".($stops[$i]+$window_size)."-".($stops[$i]+500)."|".$str;
                    push @tmp_features,$enst."|noncoding_proxintron|".$chr."|".($starts[$i+1]-500)."-".($starts[$i+1]-$window_size)."|".$str;
                    push @tmp_features,$enst."|noncoding_distintron|".$chr."|".($stops[$i]+500)."-".($starts[$i+1]-500)."|".$str;
                } else {
                    my $midpoint = int(($starts[$i+1]+$stops[$i])/2);

                    if ($starts[$i+1]-$stops[$i] > 2 * $window_size) {
                        push @tmp_features,$enst."|noncoding_proxintron|".$chr."|".($stops[$i]+$window_size)."-".($starts[$i+1]-$window_size)."|".$str;

                        if ($str eq "+") {
                            push @tmp_features,$enst."|noncoding_5ss|".$chr."|".($stops[$i])."-".($stops[$i]+$window_size)."|".$str;
                            push @tmp_features,$enst."|noncoding_3ss|".$chr."|".($starts[$i+1]-$window_size)."-".$starts[$i+1]."|".$str;
                        } elsif ($str eq "-") {
                            push @tmp_features,$enst."|noncoding_3ss|".$chr."|".($stops[$i])."-".($stops[$i]+$window_size)."|".$str;
                            push @tmp_features,$enst."|noncoding_5ss|".$chr."|".($starts[$i+1]-$window_size)."-".$starts[$i+1]."|".$str;
                        }
                    } else {
                        if ($str eq "+") {
                            push @tmp_features,$enst."|noncoding_5ss|".$chr."|".($stops[$i])."-".($midpoint)."|".$str;
                            push @tmp_features,$enst."|noncoding_3ss|".$chr."|".($midpoint)."-".$starts[$i+1]."|".$str;
                        } elsif ($str eq "-") {
                            push @tmp_features,$enst."|noncoding_3ss|".$chr."|".($stops[$i])."-".($midpoint)."|".$str;
                            push @tmp_features,$enst."|noncoding_5ss|".$chr."|".($midpoint)."-".$starts[$i+1]."|".$str;
                        }
                    }
                }
            }
        }


        for my $feature (@tmp_features) {
            my ($enst,$type,$rchr,$region,$rstr) = split(/\|/,$feature);
	    push @{$all_features{$type}},$feature;

        }
    }
    close(F);

}


sub read_gencode_gtf {

    my $file = shift;
    my $all_protein_coding_flag = shift;
#    my $file = "/projects/ps-yeolab4/genomes/hg19/gencode_v19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf";
    print STDERR "Reading in $file\n";
    open(F,$file);
    for my $line (<F>) {
        chomp($line);
        next if ($line =~ /^\#/);
        my @tmp = split(/\t/,$line);

        my $stuff = $tmp[8];
        my @stufff = split(/\;/,$stuff);
        my ($ensg_id,$gene_type,$gene_name,$enst_id,$transcript_type);

        for my $s (@stufff) {
            $s =~ s/^\s//g;
            $s =~ s/\s$//g;

            if ($s =~ /gene_id \"(.+?)\"/) {
                if ($ensg_id) {
                    print STDERR "two ensg ids? $line\n";
                }
                $ensg_id = $1;
            }
            if ($s =~ /transcript_id \"(.+?)\"/) {
                if ($enst_id) {
                    print STDERR "two enst ids? $line\n";
                }
                $enst_id = $1;
            }
            if ($s =~ /gene_type \"(.+?)\"/) {
                if ($gene_type) {
                    print STDERR "two gene types $line\n";
                }
                $gene_type = $1;
                
            }

            if ($s =~ /transcript_type \"(.+?)\"/) {
                $transcript_type = $1;
            }
            if ($s =~ /gene_name \"(.+?)\"/) {
                $gene_name = $1;
            }
        }
        next unless ($enst_id);
        if (exists $enst2ensg{$enst_id} && $ensg_id ne $enst2ensg{$enst_id}) {
            print STDERR "error two ensgs for enst $enst_id $ensg_id $enst2ensg{$enst_id}\n";
        }

        $transcript_type = "unknown" unless ($transcript_type);
        $gene_name = "unknown" unless ($gene_name);
        $gene_type = "unknown" unless ($gene_type);

        if ($all_protein_coding_flag eq "all_protein_coding") {
            $gene_type = "protein_coding";
            $transcript_type = "protein_coding";
        }

        $enst2ensg{$enst_id} = $ensg_id;
        $ensg2name{$ensg_id}{$gene_name}=1;
        $ensg2type{$ensg_id}{$gene_type}=1;
        $enst2type{$enst_id} = $transcript_type;
    }
    close(F);

}

