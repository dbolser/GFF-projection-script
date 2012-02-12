#!/usr/bin/perl -w

=head1 NAME

map_features.plx - Maps GFF features via a given 'mapping GFF' file.

=head1 SYNOPSIS

map_features.plx [options] <gff-file1> [<gff-file2>, ...]

    Options (basic):
        --help  show help
        --man   show extended help

=head1 OPTIONS

=over 8

=item --gff-mapping-file, --map

REQUIRED: The 'mapping GFF' file. See DESCRIPTION for details.

=item --feature-type, --type

OPTIONAL: The feature type to use in the 'mapping GFF' file.

=item --map_all, -a

OPTIONAL: Output all features, including those that don't get mapped.

Note, features that don't map /cleanly/ are never output, but may be
reported optionally (-v).

=head1 DESCRIPTION

The 'mapping GFF' file should have features with IDs that we want to
map /from/ and reference sequences (SEQ_IDs) that we want to map /to/.

For example, the lines of the 'mapping GFF' file may look something
like this:
chr1	src_x	contig	   1	1001	.	-	.	ID=contig1
chr1	src_x	contig	1010	3001	.	+	.	ID=contig2
chr2	src_x	contig	  20	2001	.	+	.	ID=contig3

Using this 'mapping GFF' file, we can map features /from/ contig1 or contig2
/to/ chr1, or /from/ contig3 /to/ chr2, etc.

i.e.
contig1	src_y	gene	 401	 901	3	+	1	ID=gene1
contig3	src_y	gene	  90	1409	9	-	2	ID=gene2

becomes:
chr1	src_y	gene	101	 601	3	-	1	ID=gene1
chr2	src_y	gene	109	1428	9	-	2	ID=gene2


Features on reference sequences that don't match an ID in the mapping
file, or that lay outside the region defined in the mapping file, are
dropped, or passed through unchanged, depending on the setting of
--map-all.

Features that don't map 'clenly' are dropped, and optinally reported (-v).

Finally, sub-features are 'orphaned' if their parent feature is found
not to have been mapped cleanly. To do this in a single pass, we
require parent features to preceed sub-featues in the GFF. 

=cut



use strict;

## For debugging
use Data::Dumper;

## To parse command line options, or barf
use Getopt::Long;
use Pod::Usage;

## To parse GFF3
use Bio::GFF3::LowLevel qw/ gff3_parse_feature gff3_format_feature /;

## BioPerl 'location' features
use Bio::Location::Simple;

## Our Moose powered wrapper to Bio::Coordinate::Collection
use GFFCoordinateMapper;



## COMMAND LINE OPTIONS

my $man = 0;
my $help = 0;
my $debug = 0;
my $verbose = 0;

my $gff_mapping_file;
my $feature_type;

my $map_all = 0;

GetOptions( 'help|?'   => \$help,
            'man'      => \$man,
            'debug+'   => \$debug,
            'verbose+' => \$verbose,

            'gff-mapping-file|map=s' => \$gff_mapping_file,
            'feature-type|type=s'    => \$feature_type,

            'map-all-features|a+' => \$map_all,
          )
  or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -verbose => 2 ) if $man;

pod2usage( "\nplease pass a 'mapping' GFF file (-g)\n" )
  unless -s $gff_mapping_file;

pod2usage( "\nPass a GFF file (or files) for me to map!\n" )
  unless @ARGV;

warn "debug\t: $debug\nverbose\t: $verbose\n"
  if $verbose > 0;





## GO TIME

warn "\nprocessing mapping GFF\n";

## Initalise the 'mapping' object
my $mapper;

if($feature_type){
  $mapper = GFFCoordinateMapper->
    new( file => $gff_mapping_file,
	 type => $feature_type
       );
}
else{
  $mapper = GFFCoordinateMapper->
    new( file => $gff_mapping_file );
}

die "found zero components to map over in $gff_mapping_file\n\n"
  unless $mapper->components > 0;

warn "the mapper has ",
  scalar $mapper->components, " components\n\n";



warn "processing GFF\n";

my ($features_mapped,
    %features_mapped,
    $failed_to_map,
    %failed_to_map);

## NOTE, below we expect all sub-features to be preceeded by their
## parent features in the GFF, because we 'orphan' features who's
## parents didn't get mapped.

while(<>){
    next if /^#/;
    next if /^\s*$/;
    
    ## Parse the GFF
    my $feature =
      gff3_parse_feature( $_ );
    
    ## Debugging
    if($debug > 0){
        next unless
          $feature->{attributes}{ID}[0] =~ /PGSC0003DMG400001256/;
    }
    
    ## Debugging
    warn Dumper $feature
      if $debug > 2;
    
    ## Build a location object from the feature
    my $feature_location = Bio::Location::Simple->
      new( -seq_id => $feature->{seq_id},
           -start  => $feature->{start},
           -end    => $feature->{end},
           -strand => $feature->{strand},
         );
    
    warn Dumper $feature_location
      if $debug > 2;
        
    
    
    ##
    ## Try to map the feature location onto the new coordinates
    ##
    
    my $new_feature_location =
      $mapper->map( $feature_location );
    
    warn Dumper $new_feature_location
      if $debug > 2;
    
    
    
    ## After mapping, there are three possibilities:
    
    ## 1) The feature maps outside the coordinate space defined by the
    ##    mapper and should be ignored or passed through unchanged.
    
    ## 2) The feature maps cleanly, and its coordinates should be
    ##    updated to their new position.
    
    ## 3) The feature may span multiple regions, and should now be
    ##    dropped.
    
    ## The mapping result object allows us to conveniently differentiate
    ## between these three cases...
    
    my $num_gaps  = $new_feature_location->each_gap;
    my $num_match = $new_feature_location->each_match;
    
    
    
    if(0){} # I hate syntax
    
    ## Case 1, nothing much to do here
    elsif($num_gaps == 1 && $num_match == 0){
        next unless $map_all;
    }
    
    ## Case 2, adjust the feature coordinates
    elsif($num_gaps == 0 && $num_match == 1){
        $feature->{seq_id} = $new_feature_location->seq_id;
        $feature->{start}  = $new_feature_location->start;
        $feature->{end}    = $new_feature_location->end;
        $feature->{strand} = $new_feature_location->strand;
    }
    
    ## Case 3, something else happened
    else{
        warn "spanning feature! : ",
          $feature->{attributes}{ID}[0], "\n"
            if $verbose > 0;
        
        warn Dumper $feature
          if $debug > 1;
        warn Dumper $new_feature_location
          if $debug > 1;
        
        ## log our 'failure'...
        $failed_to_map++;
        push @{$failed_to_map{$feature->{seq_id}}}, $feature;
        
        ## and move on without printing
        next;# unless $really_map_all;
    }
    
    
    
    ## Made it!
    $features_mapped++;
    $features_mapped{$feature->{attributes}{ID}[0]}++;
    
    
    
    ## If this feature has a parent feature, check that the parent
    ## feature made it through, else, orphan the feature (so
    ## cruel!). Note, although this choice doesn't make much sense for
    ## gene sub-features, it's important for 'clone-end' features. It
    ## could be made optional, or conditional on feature type.
    
    if($feature->{attributes}{Parent}){
        my @parents = @{$feature->{attributes}{Parent}};
        my @parents_mapped = grep $features_mapped{$_}, @parents;
        if(@parents_mapped < @parents){
            warn "removing unmapped parent from ",
              $feature->{attributes}{ID}[0], "\n"
                if $verbose > 0;
            if(@parents_mapped){
                $feature->{attributes}{Parent} = \@parents_mapped;
            }
            else{
                delete $feature->{attributes}{Parent};
            }
        }
    }
    
    
    
    ## Convert from the BioPerl standard to the GFF standard... Sigh...
    if($feature->{strand}){
        $feature->{strand} = '+' if $feature->{strand} eq  '1';
        $feature->{strand} = '-' if $feature->{strand} eq '-1';
        $feature->{strand} = '?' if $feature->{strand} eq  '0';
    }
    
    ## Write it
    print gff3_format_feature( $feature );
}

warn "we mapped ", $features_mapped || 0, " features\n";
warn "there were ", $failed_to_map || 0, " problem features\n\n";



warn "error reporting : \n\n";

for (sort keys %failed_to_map){
    warn 'failed to map ',
      scalar @{$failed_to_map{$_}}, " features on $_\n";
    
    for(@{$failed_to_map{$_}}){
        warn "\t", $_->{attributes}{ID}[0], "\t",
          ($_->{attributes}{Name}[0] || ''), "\n"
            if $verbose > 0;
        warn Dumper $_
          if $debug > 1;
    }
}

warn "OK\n";
