=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>
    
=cut

=head1 NAME

 NearestExon

=head1 SYNOPSIS

 mv NearestExon.pm ~/.vep/Plugins
 ./vep -i variations.vcf --cache --plugin NearestExon

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 finds the nearest exon junction boundary(ies) to a variant. More than one boundary
 may be reported if the boundaries are equidistant.
 
 This plugin does not run in offline mode.

 Various parameters can be altered by passing them to the plugin command:

 - limit     : limit the number of exons returned (default: 1)
 - range     : initial search range in bp (default: 1000)
 - max_range : maximum search range in bp (default: 10000)

 Parameters are passed e.g.:

 --plugin NearestExon,limit=3,max_range=50000

=cut

package NearestExon;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my %CONFIG = (
  limit => 1,
  range => 1000,
  max_range => 10000,
);

my $char_sep = "|";

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $params = $self->params;

  # get output format
  $char_sep = "+" if ($self->{config}->{output_format} eq 'vcf');

  foreach my $param(@$params) {
    my ($key, $val) = split('=', $param);
    die("ERROR: Failed to parse parameter $param\n") unless defined($key) && defined($val);
    $CONFIG{$key} = $val;
  }

  die("ERROR: This plugin does not work in --offline mode\n") if $self->{config}->{offline};

  return $self;
}

sub feature_types {
  return ['Feature','Intergenic'];
}

sub variant_feature_types {
  return ['BaseVariationFeature'];
}

sub get_header_info {
  my $header = 'Nearest Exon. Format:';
  $header .= join($char_sep, qw(ExonID distance start/end) );

  return {
    NearestExon => $header,
  }
}

sub run {
  my ($self, $vfoa) = @_;

  my $vf = $vfoa->base_variation_feature;
  my $loc_string = sprintf("%s:%i-%i", $vf->{chr} || $vf->seq_region_name, $vf->{start}, $vf->{end});

  if(!exists($self->{_cache}) || !exists($self->{_cache}->{$loc_string})) {
    $self->{config}->{ea} = $self->{config}->{reg}->get_adaptor($self->{config}->{species}, $self->{config}->{core_type}, 'Exon');
    $self->{ea} ||= $self->{config}->{ea};
    die("ERROR: Could not get exon adaptor;\n") unless $self->{ea};

    my %opts = map {'-'.$_ => $CONFIG{$_}} keys %CONFIG;
    $opts{-feature} = $vf;

    my @result = map {$_->[0]} @{
      $self->{ea}->fetch_all_by_outward_search(%opts)
    };

    my %dists;
    my $min = $CONFIG{max_range};
    foreach my $exon (@result){
      my $startD = abs ($vf->start - $exon->seq_region_start);
      my $endD = abs ($vf->start - $exon->seq_region_end);
      if ($startD < $endD){
        $dists{$exon->stable_id}{$startD} = 'start';
        $min = $startD if $min > $startD;
      } elsif ($startD > $endD){
        $dists{$exon->stable_id}{$endD} = 'end';
        $min = $endD if $min > $endD;
      } else {
        $dists{$exon->stable_id}{$startD} = "start_end";
        $min = $startD if $min > $startD;
      }
    }

    my @finalRes;
    foreach my $exon (keys %dists){
      if (exists $dists{$exon}{$min}) {
        push(@finalRes, $exon.$char_sep.$min.$char_sep.$dists{$exon}{$min})
      }
    }

    $self->{_cache}->{$loc_string} = scalar @finalRes ? join(',', @finalRes) : undef;
  }

  return $self->{_cache}->{$loc_string} ? { NearestExon => $self->{_cache}->{$loc_string} } : {};
}

1;

