package GFFCoordinateMapper;

use Data::Dumper;

use Moose;
use MooseX::FileAttribute;

use Bio::Location::Simple;
use Bio::Coordinate::Pair;
use Bio::Coordinate::Collection;

use Bio::GFF3::LowLevel qw/ gff3_parse_feature /;


has 'mapper' =>
  (
   is => 'ro',
   isa => 'Bio::Coordinate::Collection',

   ## Delegate
   handles => [ qw( map add_mapper swap ) ],

   ## Build it
   default => sub{ Bio::Coordinate::Collection->new() },
  );

## Just so that we can pass a file at create time...
has_file 'gff_file' =>
  (
   must_exist => 1,
   init_arg => 'file',
   trigger => \&load_gff,
  );

has 'feature_type' =>
  (
   is => 'ro',
   isa => 'Str',
   init_arg => 'type',
  );



=head1 NAME

GFFCoordinateMapper - Perl module for building a
                      Bio::Coordinate::Collection from a GFF file.

=head1 SYNOPSIS

  use GFFCoordinateMapper;

  # Create a GFF coordinate mapper, and pass it a GFF file to build
  # mappings from
  $feature_mapper = GFFCoordinateMapper->new();
  $feature_mapper->load_gff( $gff_file );

  # As above, but in a single step
  $feature_mapper = GFFCoordinateMapper->
    new( file => $mapping_gff );

  # As above, but define a specifc feature to use for the mapping.
  # The default is to use all features.
  $feature_mapper = GFFCoordinateMapper->
    new( file => $mapping_gff, type => 'type_of_feature' );

=head1 DESCRIPTION

This module builds a Bio::Coordinate::Collection (a mapper) from the
given GFF file. A mapper is a collection of Bio::Coordinate::Pair
objects. The resulting object is suitable for 'feature mapping'
through the ID to SEQ_ID relationships defined in the GFF.

i.e. If the GFF describes a contig (ID=cx) on a reference chromosome
(SEQ_ID=chr1), you can use the resulting object to map features with
coordinates on cx (SEQ_ID=cx) to their propper position on chr1.

See Bio::Coordinate::Collection.

=head2 Methods

=over 4

=item * load_gff

Takes a GFF file as an argument and uses it to build the
Bio::Coordinate::Collection. If feature_type is set, only those
features with that type will be used.

Returns the (number of) Bio::Coordinate::Pair 'components' in the
mapper.

=item * components

A convenience method.

Returns the (number of) Bio::Coordinate::Pair 'components' in the
mapper.

=back

=head2 Methods from Bio::Coordinate::Collection

=over 4

=item * map

Takes a SeqFeature and tries to map it through the mapper. It uses the
SEQ_ID, START and END of the feature to match an appropriate
Bio::Coordinate::Pair and perform the mapping.

Returns a feature location object (IIRC)

    After mapping, there are three possibilities:

    1) The feature maps outside the coordinate space defined by the
       mapper. The feature location is a single gap.

    2) The feature maps cleanly. The feature location is a single
       match.

    3) The feature may span multiple regions. The feature location is
       a mix of more than one gap and match.

TODO: Return the mapped feature, not the resulting feature location!

=item * add_mapper

Adds another Bio::Coordinate::Pair to the Bio::Coordinate::Collection

=item * swap

Swaps the direction of mapping, such that you can go from chromosomes
to contigs or vice verse.

=back

=head1 AUTHOR

Dan B (dan.bolser@gmail.com)

=cut



sub load_gff {
  my $self = shift;
  my $file = shift;

  confess "pass a GFF file plz\n"
    unless -s $file;

  open C, '<', $file
    or die "failed to open file '$file' : $!\n";

  while(<C>){
    ## Ignore comments or blank lines
    next if /^#/;
    next if /^\s*$/;

    ## Parse the GFF
    my $feature =
      gff3_parse_feature( $_ );

    ## If feature_type is set, respect its wishes
    if ($self->feature_type){
      next unless $feature->{type} eq $self->feature_type;
    }

    ## Get the 'extent' of the feature
    my $feature_length = $feature->{end} - $feature->{start} + 1;

    ## Create a Bio::Coordinate::Pair (map) to store the mapping
    ## between the feature and its reference sequence

    my $seq_id_loc = Bio::Location::Simple->
      new( -seq_id => $feature->{seq_id},
	   -start  => $feature->{start},
	   -end    => $feature->{end},
	   -strand => +1,
	 );

    my $featur_loc = Bio::Location::Simple->
      new( -seq_id => $feature->{attributes}{ID}[0],
	   -start  => 1,
	   -end    => $feature_length,
	   -strand => $feature->{strand},
	 );

    my $map = Bio::Coordinate::Pair->
      new( -in  => $featur_loc,
	   -out => $seq_id_loc,
	 );

    $self->add_mapper( $map );
  }

  return $self->components;
}

sub components {
  my $self = shift;
  
  my @components;
  
  ## $self->agp is a Bio::Coordinate::Colleciton, composed of several
  ## 'mappers'. Each 'mapper' is a Bio::Coordinate::Pair.
  push @components,
  $_->in->seq_id for $self->mapper->mappers;
  
  return @components;
}

1;
