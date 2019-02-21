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

 NearestGene

=head1 SYNOPSIS

 mv NearestExon.pm ~/.vep/Plugins
 ./vep -i variations.vcf --cache --plugin NearestExonJB

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 finds the nearest exon junction boundary(ies) to a variant. More than one boundary
 may be reported if the boundaries are equidistant.

 Various parameters can be altered by passing them to the plugin command:

 - max_range : maximum search range in bp (default: 10000)

 Parameters are passed e.g.:

 --plugin NearestExon,limit=3,max_range=50000

=cut

package NearestExonJB;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use Data::Dumper;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my $char_sep = "|";

my %CONFIG = (
  max_range => 10000,
);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $params = $self->params;

  # get output format
  my $out_vcf  = 1 if ($self->{config}->{output_format} eq "vcf");
  my $out_json = 1 if ($self->{config}->{output_format} eq "json");
  my $out_txt = 0 if ($out_vcf || $out_json);
  $char_sep = "+" if $out_vcf;

  foreach my $param(@$params) {
    my ($key, $val) = split('=', $param);
    die("ERROR: Failed to parse parameter $param\n") unless defined($key) && defined($val);
    $CONFIG{$key} = $val;
  }

  return $self;
}

sub feature_types {
  return ['Feature','Intergenic'];
}

#TODO: do I need this?
sub variant_feature_types {
  return ['BaseVariationFeature'];
}


sub get_header_info {
  my $header = 'Nearest Exon Junction Boundary. Format:';
  $header .= join($char_sep, qw(ExonID distance start/end) );

  return {
    NearestExonJB => $header,
  }
}

sub run {
  my ($self, $vfoa) = @_;

  my $vf = $vfoa->base_variation_feature;
  my $trv = $vfoa->base_transcript_variation;
  my $loc_string = sprintf("%s:%i-%i", $vf->{chr} || $vf->seq_region_name, $vf->{start}, $vf->{end});

  if(!exists($self->{_cache}) || !exists($self->{_cache}->{$trv->transcript_stable_id})) {
    $self->{config}->{ea} = $self->{config}->{reg}->get_adaptor($self->{config}->{species}, $self->{config}->{core_type}, 'Exon');
    $self->{ea} ||= $self->{config}->{ea};
    die("ERROR: Could not get exon adaptor; this plugin does not work in --offline mode\n") unless $self->{ea};

    my @exons = @{$trv->transcript->get_all_Exons};
    my %dists;
    my $min = $CONFIG{max_range};

    if (defined $trv->exon_number){
      my @tmp = split('/',$trv->exon_number);
      my $exon = $exons[$tmp[0]-1];
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

    $self->{_cache}->{$trv->transcript_stable_id} = scalar @finalRes ? join(",", @finalRes) : undef;
  }
  return $self->{_cache}->{$trv->transcript_stable_id} ? { NearestExonJB => $self->{_cache}->{$trv->transcript_stable_id} } : {};
}

1;

