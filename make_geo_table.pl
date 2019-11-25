#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use LWP::Simple;
use HTML::TableExtract;
use HTML::TokeParser;

#********************************************
#author: laurieCannon 
#program: make_geo_table.pl
#description: make a geo table!
#*******************************************
#Summary/Program Construction:
#
#input file: list of accessionNum geoUrl; tab delimited 
#label tab delimited table/column headers, \n
#read input file as FH, by line
#write line 1 of input file: accessionNum, \t, geoUrl
#wget geoUrl > output to wget_html directory
#Slurp html output file as string
#parse output html file for geo table information
#use regex of html <td> (table data) tags to capture infomation 
#delete html output file when all geo table info captured from it
#read/write input file line 2: accessionNUm, \t, geoUrl
#...blahblah until end of inputFile
#*********************************************

my $inputFile; #userdefined accessionNum and geoUrl file, tab delimited
my $outputFile; #userdefined geo table output file

GetOptions(
            'i:s' => \$inputFile,
            'o:s' => \$outputFile,
            ) or die $!;
            

open (OUTPUT, ">", $outputFile) or die $!;
print OUTPUT "Accession","\t","URL","\t","Release Date","\t","Title","\t","Organism","\t","Experiment Type","\t","Summary","\t","Overall Design","\t","Assays","\t","Processed","\t","RawData","\t","PresentInAtlas","\t","ArrayExpressURL","\n";

#retrieve geo content
my $url; #url address from input file at each line
my $html; #retrieve html data from current FH line

my $htmlFile = "html.txt"; #concat all of the html 
my @line;

open ( FH, $inputFile) or die "no file: $!";
    
    my $next;
    while ($next = <FH>){
   # chomp($next); 
 open (HTML_OUTPUT, ">", $htmlFile) or die $!;
 print HTML_OUTPUT "file";
    @line = split('\t', $next);
  
    print $line[0], "\t",$line[1], "\n";
    $url = $line[1];
    $html = get($url) or die "unable to get URL page: $!";
    print HTML_OUTPUT "$html";
     
my $token;
my $nextToken;
my $textBody;

#my $platforms =~ (m/"Platforms".{}/g); #Platforms tag varies across html: i.e. "Platforms (2)", "Platforms (3)"

#my $samples = ($html =~ m/^("<td>Samples")/); #see $platforms for reasoning
my $x=0; #two sample td tags per geo page: want second one



$line[1] =~ s/\r\n//g; #remove carriage return/space after url

print OUTPUT "$line[0]\t$line[1]\t";

my $sampleText;
my $sampleTag;
my @row;
push(@row, $line[0]);
push (@row, $line[1]);

my  $p = HTML::TokeParser->new(shift || "html.txt");

  while ($token = $p->get_tag("td")) {
       my $text = $p->get_trimmed_text("/td");
       if($text =~ m/^(Samples)(.)/g){
       print "MATCHED\n";
       
            print "sampleText 1\n"; 
                if($x==1){ #get the sample information
                    
                    print "$sampleText 2\n"; 
                    while( ($sampleTag = $p -> get_tag("td")) ne "SRA" ){
                         
                            $sampleText = $p -> get_trimmed_text("/td");
                            print "$sampleText\n";   
                    }
                }   #end of $x==1 get sample info
                else{ 
                    $x++;
                    }#the next $sample capture will be the $sample tag we want, when $x==1
        }#end of if $text = $samples
            
  }#end while $token
        
            $nextToken = $p->get_tag("td");
            $textBody = $p->get_trimmed_text("/td");
            print OUTPUT "$textBody\t";
            push(@row, $textBody);
         
   
    print OUTPUT "\n";
  
   

   
   








