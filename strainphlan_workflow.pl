#!/usr/bin/perl
use strict;
use FindBin qw($Bin $Script);
use lib $Bin;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename;

my $cwd = abs_path;
my $cmd = basename $0;

sub usage {
        print <<USAGE;
usage:
        $cmd
pattern
        -id|sample_fq:
        -phe|metadata:
        -ref|bac_ref:
        -F|remarker:
        -pf|prefix:
        -cmd|command:  T(true,defalut) or F(false), if False, don't to excute the step1 and step2;
USAGE
};

my $pattern = $ARGV[0];
print &usage && exit if(!defined $pattern);


my ($sample_fq, $metadata, $bac_ref, $remarker, $prefix);

GetOptions(
                "id:s" => \$sample_fq,
                "phe:s" => \$metadata,
                "ref:s" => \$bac_ref,
                "F:s" => \$remarker,
                "pf:s" => \$prefix,
        );



my $workDir = $prefix."_strain";
my $currentDir = abs_path(".");
my @sample_name; # sample id

system "mkdir -p $workDir" unless(-d $workDir);
$workDir = abs_path($workDir);


# step1 Run MetaPhlAn2
open SH,">$workDir/01.metaphlan2.sh" || die $!;

my $step1_res = $workDir."/metaphlan";
system "mkdir -p $step1_res" unless(-d $step1_res);


open I, "$sample_fq" or die $!;
while(<I>){
        chomp;
        my $name = (split /\//, $_)[-1];
        my $sample ;
        if($name=~/fq.gz/){
                $name=~/(\S+).fq.gz/;
                $sample = $1;
        }else{
                $name=~/(\S+).fq/       ;
                $sample = $1;
        }
        push(@sample_name, $sample);
        print SH "metaphlan2.py $_ $step1_res/$sample\_profile.txt --bowtie2out $step1_res/$sample\_bowtie2.txt --samout $step1_res/$sample.sam.bz2 --input_type fastq\n";
}

close SH;

# step2 sam2marker
# generate a marker file for each sample. The marker files contain the consensus of unique marker genes for each species found in the sample, which will be used for SNP profiling. Run the commands below (with the caveat above!) to generate the marker files in your current working directory.


open SH1, ">$workDir/02.sam2marker.sh" || die $!;
my $step2_consensus_markers = $workDir."/consensus_mark";
system "mkdir -p $step2_consensus_markers" unless(-d $step2_consensus_markers);

foreach my $id (@sample_name){
        print SH1 "source activate strainphlan \nsample2markers.py --ifn_samples $step1_res/$id.sam.bz2 --input_type sam --output_dir $step2_consensus_markers --nprocs 8 &>$step2_consensus_markers/log.txt\n";
        }

close SH1;

# step3 extract_db_marker
# Identify clades detected in the samples and build reference databases

open SH2, ">$workDir/03.extract_db_marker.sh"|| die $!;
my $db_marker = "/hwfssz1/ST_META/share/database/strainphlan_db_markers";
my $mpa_db = "/ldfssz1/ST_META/share/flow/biobakery-metaphlan2-d8ab9ca4244c";

print SH2 "source activate strainphlan \nstrainphlan.py --ifn_samples $step2_consensus_markers/*markers --output_dir $workDir --print_clades_only --nprocs_main 8 >clades.txt\n";

#if(-e "$workDir/clades.txt"){
open I2,"$currentDir/clades.txt" or die $!;
my @bac_list;
while(<I2>){
        chomp;
        push(@bac_list, $_);
}

if(grep {$bac_ref == $_} @bac_list){
                print SH2 "extract_markers.py --mpa_pkl $mpa_db/db_v20/mpa_v20_m200.pkl --ifn_markers $db_marker/all_markers.fasta --clade $bac_ref --ofn_markers db_markers/$bac_ref.fasta ";
        }else{
                print "this bac is not in the clade\n";
        #       die $!;
        }
#}else{
#       next;
#}
# step4 generate trees from alignments
# step5 plot,Visualization with ggtree

# Run StrainPhlAn to generate alignments and then trees, providing the sample marker files and
# the clade reference marker file generated in prior steps

# add the ref, here use the panphlan database, if the bac_ref not in the db ,should creat manully

my $ggtree = "/hwfssz1/ST_META/share/User/renhuahui/breadcrumbs/breadcrumbs/scripts";


#open SH3, ">$workDir/04.build_metadata.sh" or die $!;
open SH4, ">$workDir/04.build_tree.sh" or die $!;
open SH5, ">$workDir/05.drawing.sh" or die $!;

my $output = $workDir."/output";
my $metadir = $workDir."/metadata";
my $img = $workDir."/img";
system "mkdir -p $output" unless(-e $output);
system "mkdir -p $metadir" unless(-e $metadir);
system "mkdir -p $img" unless(-e $img);

my %bac_ref;
open I3, "/hwfssz1/ST_META/share/database/panphlan_fnn/Bac_ffn2.list" or die $!;
while(<I3>){
        chomp;
        my @tmp_str = split /\t/, $_;
        $bac_ref{$tmp_str[0]} = $tmp_str[1];

}


# here should add the ref to  the metadata


foreach my $i (@bac_list){
        my @str = split /\s+/,$i;
        my $len = @str;
        my $bac;
        #print "$len";
        if($len == 2){
                my $str1 = $str[1];
        #       print "$str1";
                if($str1=~/\|/){
                                my @str2 = split /\|/,$str1;
                                my $tmp = $str2[-1];
                                $tmp =~ /(\S+)\)/;
                                $bac = $1;
                                #print "$str1\n"        ;
                        }else{
                                $str1=~/\((\S+)\)/;
                                $bac = $1;
                                #print "$bac\n";
                        }
        }else{
                $bac = $i;
                }
        #print "$bac\n";
        if(exists $bac_ref{$bac}){
        print "$bac_ref{$bac}\n";

        # generate metadata.txt
        open I4,">$metadir/metadata.$bac.txt";
        open I5,"$metadata" or die $!;
        while(<I5>){
                print I4 "$_";
                }
        my $tmp_name = (split /\//,$bac_ref{$bac})[-1];
        $tmp_name=~/(\S+).ffn/;
        print I4 "$1\tRefGenome\n";

        print SH4 "source activate strainphlan \nstrainphlan.py --mpa_pkl $mpa_db/db_v20/mpa_v20_m200.pkl --ifn_samples $step2_consensus_markers/*markers --ifn_markers $db_marker/$bac.markers.fasta --ifn_ref_genomes $bac_ref{$bac} --output_dir $output --nprocs_main 10 --clades $str[0] &> output/$bac._log_full.txt\n";

        print SH4 "$mpa_db/strainphlan_src/add_metadata_tree.py --ifn_trees $output/RAxML_bestTree.$bac.tree --ifn_metadatas $metadir/metadata.$bac.txt --metadatas subjectID\n";

        print SH5 "$ggtree/strainphlan_ggtree.R $output/RAxML_bestTree.$bac.tree metadata.txt $output/$bac.fasta $img/$bac._tree_1.png $img/$bac._tree_2.png\n" ;

        close I4;
        close I5;
        }else{
                next;
        }
}
# step6  snp
