#!/usr/bin/env perl
use strict;
use warnings;

my $version = '
Version: 0.1.1 (2026-01-29)
Author: Adam Yongxin Ye @ BCH & ChatGPT
';

my $usage = "Usage: perl $0 <input.tlx> <ref.fa>
  [CAC_slop_bp (default:15)]
  [bait_hint (default:2000)]
  [major_bait_CAC_reuse_max_dist (default:200)]
  [ideal_vs_cryptic_dist_thres (default:15)]
  [max_PelementLen (default:4)]
  [CAC_search_window (default:-3,200)]
  [debug_Qname (default:none)]


Input:
  <input.tlx>    TLX format with header line.
    Must contain columns: Qname, Seq,
      Rname, Rstart, Rend, Strand, Qstart,
      B_Rname, B_Rstart, B_Rend, B_Strand, B_Qend
  <ref.fa>       Reference FASTA.
    Only TLX records with Rname and B_Rname both exist in ref.fa will be output.

Options:
  CAC_slop_bp (default:15)
    Number of bp extracted around CAC when reporting P_CAC_ref_seq / B_CAC_ref_seq.
    Does not affect CAC searching.
  bait_hint (default:2000)
    Either an integer N (to use first N TLX records to infer majority bait CAC),
    or chr:pos for majority bait CAC, where the sign of pos indicates strand.
  major_bait_CAC_reuse_max_dist (default:200)
    Use majority bait CAC if bait junction is within this distance (same chr and strand).
  ideal_vs_cryptic_dist_thres (default:15)
    Threshold to choose ideal vs cryptic motif when both are found near a junction.
    if ideal_delLen-cryptic_delLen > ideal_vs_cryptic_dist_thres, choose cryptic; otherwise, choose ideal
  max_PelementLen (default:4)
    Maximum palindromic element length checked for P_PelementLen and B_PelementLen.
  CAC_search_window (default:-3,200)
    CAC search window defined as intrude_bp,search_bp.
    Allows prey alignment to intrude CAC motif by abs(intrude_bp) bp
    and searches search_bp bp outward from junction.
  debug_Qname (default:none)
    If provided, print extra debug information for the matching Qname to STDERR.

Output:
  STDOUT. TLX format with additional columns appended:
    midLen = Qstart-B_Qend-1    length of mid (read sequence between bait and prey)
      can represent insertion (+) or microhomology (-)
    B_CAC_coord (chr:pos with sign for strand)    determined first C of CAC nearby bait junction
      CAC nearby prey and bait junctions is determined on strand-specific rules
    B_CAC_ref_seq    reference sequence at bait CAC (+-CAC_slop_bp), CAC... uppercase, upsteram lowercase
    B_delLen    deletion length (distance) between bait junction and the first C of nearby CAC
      can be negative if greedy alignment (reported by original HTGTS pipeline) happens to intrude into CAC
    B_delLen_adj = max(0, B_delLen)    bait deletion length adjusted for negative B_delLen
    P_CAC_coord (chr:pos with sign for strand)    determined first C of CAC nearby prey junction
    P_CAC_ref_seq    reference sequence at prey CAC (+-CAC_slop_bp), CAC... uppercase, upsteram lowercase
    P_delLen    deletion length (distance) between prey junction and the first C of nearby CAC
    P_delLen_adj = max(0, P_delLen)   prey deletion length adjusted for negative P_delLen
    B_Qend_adj    bait junction on read adjusted for negative B_delLen (intruding to CAC)
    B_Junction_adj    bait junction on ref adjusted for negative B_delLen (intruding to CAC)
    Qstart_adj    prey junction on read adjusted for negative P_delLen (intruding to CAC)
    Junction_adj    prey junction on ref adjusted for negative P_delLen (intruding to CAC)
      when delLen < 0, adjust to force delLen to 0 and shift junction coordinate to the base just before CAC
    midLen_adj = Qstart_adj-B_Qend_adj-1    length of mid adjusted for negative P_delLen and B_delLen
    insLen = max(0, midLen_adj)   length of insertion
    P_PelementLen    length of palindromic insertion at prey junction (greedy up to max_PelementLen)
    B_PelementLen    length of palindromic insertion at bait junction (greedy up to max_PelementLen)
      P-insertions are detected only when delLen_adj = 0 and are capped by max_PelementLen.
    nonPinsLen = max(0, midLen_adj-P_PelementLen-B_PelementLen)    length of N-insertion
    MHlen = max(0, -(midLen_adj-P_PelementLen-B_PelementLen))    length of microhomology

Notes:
  TLX coordinates are 1-based and end-included.
  Strand values (1, -1) are interpreted by numeric sign (>0 or <0).
  Missing or undefined values are outputted as -.
  Qstart = the start position of prey on read (reported by original HTGTS pipeline, after its greedy alignment to reference)
  B_Qend = the end position of bait on read (reported by original HTGTS pipeline, after its greedy alignment to reference)
  Junction = the position of prey junction on reference = Rstart if Strand > 0; or Rend if Strand < 0
  B_Junction = the position of bait junction on reference = B_Rend if B_Strand > 0; or B_Rstart if BStrand < 0
  When delLen < 0, junctions are adjusted to be adjacent to CAC and delLen_adj becomes 0.
    define intrudeRSSlen = -raw_delLen, then adjust read coords:
      B_Qend_adj = B_Qend - B_intrudeRSSlen
      Qstart_adj = Qstart + P_intrudeRSSlen
  If P_PinsLen + B_PinsLen exceeds midLen_adj, PinsLens are reduced until sum to midLen_adj,
    yielding P_PinsLen_adj, B_PinsLen_adj, and NinsLen = 0.
".$version;

# -------------------------------
# Args (positional; no Getopt)
# -------------------------------
if(@ARGV < 2){
  die $usage;
}
my $in_tlx = shift(@ARGV);
my $ref_fa = shift(@ARGV);

my $CAC_slop_bp = 15;
$CAC_slop_bp = shift(@ARGV) if (@ARGV > 0);
die "Error: CAC_slop_bp must be a positive integer.\n" if (!defined $CAC_slop_bp || $CAC_slop_bp !~ /^\d+$/ || $CAC_slop_bp < 1);

my $bait_hint = undef;
$bait_hint = shift(@ARGV) if (@ARGV > 0);

my $maj_reuse_max_dist = 200;
$maj_reuse_max_dist = shift(@ARGV) if (@ARGV > 0);
die "Error: major_bait_CAC_reuse_max_dist must be a non-negative integer.\n" if (!defined $maj_reuse_max_dist || $maj_reuse_max_dist !~ /^\d+$/);

my $ideal_vs_cryptic_dist_thres = 15;
$ideal_vs_cryptic_dist_thres = shift(@ARGV) if (@ARGV > 0);
die "Error: ideal_vs_cryptic_dist_thres must be a non-negative integer.\n" if (!defined $ideal_vs_cryptic_dist_thres || $ideal_vs_cryptic_dist_thres !~ /^\d+$/);

my $maxPinsLen = 4;
$maxPinsLen = shift(@ARGV) if (@ARGV > 0);
die "Error: max_PinsLen must be a non-negative integer.\n" if (!defined $maxPinsLen || $maxPinsLen !~ /^\d+$/);

my $CAC_search_window = "-3,200";
$CAC_search_window = shift(@ARGV) if (@ARGV > 0);

my $debug_Qname = undef;
$debug_Qname = shift(@ARGV) if (@ARGV > 0);
$debug_Qname = undef if (defined $debug_Qname && $debug_Qname eq "none");

# Parse CAC_search_window: intrude_bp,search_bp
my ($maxIntrudeIntoAlnBp, $search_win) = (3, 200);
if (defined $CAC_search_window && $CAC_search_window ne "") {
  if ($CAC_search_window =~ /^\s*(-?\d+)\s*,\s*(\d+)\s*$/) {
    $maxIntrudeIntoAlnBp = abs(int($1));
    $search_win          = int($2);
  } else {
    die "Error: CAC_search_window must be intrude_bp,search_bp (e.g. -3,200).\n";
  }
}

# Motifs on reference + strand
my $ideal_heptamer_minus_on_plus  = "CACTGTG";
my $cryptic_trimer_minus_on_plus  = "GTG";
my $ideal_heptamer_plus_on_plus   = "CACAGTG";
my $cryptic_trimer_plus_on_plus   = "CAC";

# Cache: key => [cac_strand, cpos_abs, raw_delLen]
my %cac_cache;

# -------------------------------
# Load reference fasta
# -------------------------------
my (%ref_seq, %ref_len);

open(my $fh_fa, "<", $ref_fa) or die "Error: cannot open ref.fa: $ref_fa\n";
my ($cur, $seq) = ("", "");
while (my $line = <$fh_fa>) {
  $line =~ s/[\r\n]+$//;
  if ($line =~ /^>(\S+)/) {
    if ($cur ne "") {
      $seq =~ s/\s+//g;
      $seq = uc($seq);
      $ref_seq{$cur} = $seq;
      $ref_len{$cur} = length($seq);
    }
    $cur = $1;
    $seq = "";
  } else {
    $seq .= $line;
  }
}
if ($cur ne "") {
  $seq =~ s/\s+//g;
  $seq = uc($seq);
  $ref_seq{$cur} = $seq;
  $ref_len{$cur} = length($seq);
}
close($fh_fa);
die "Error: empty reference loaded from $ref_fa\n" if (!%ref_seq);

# -------------------------------
# Helpers
# -------------------------------
sub _min { return $_[0] < $_[1] ? $_[0] : $_[1]; }
sub _max { return $_[0] > $_[1] ? $_[0] : $_[1]; }

sub outv {
  my ($v) = @_;
  return defined $v ? $v : "-";
}

sub revcomp {
  my ($s) = @_;
  $s = "" if (!defined $s);
  $s = uc($s);
  $s =~ tr/ACGTN/TGCAN/;
  return scalar reverse($s);
}

sub parse_coord {
  my ($s) = @_;
  return (undef, undef) if (!defined $s);
  if ($s =~ /^(\S+):([+-]?\d+)$/) {
    return ($1, int($2));
  }
  return (undef, undef);
}

sub format_coord {
  my ($chr, $cac_strand, $cpos_abs) = @_;
  return undef if (!defined $chr || !defined $cac_strand || !defined $cpos_abs);
  my $signed = ($cac_strand eq "-") ? -int($cpos_abs) : int($cpos_abs);
  return $chr . ":" . $signed;
}

sub strand_sign {
  my ($s) = @_;
  return 0 if (!defined $s);
  $s =~ s/^\s+|\s+$//g;
  return 0 if ($s eq "");
  my $v = int($s);
  return 0 if ($v == 0);
  return ($v > 0) ? 1 : -1;
}

sub build_cac_seq {
  my ($chr, $cpos_abs, $cac_strand, $slop) = @_;
  return undef if (!defined $chr || !defined $cpos_abs || $cpos_abs !~ /^\d+$/);
  return undef if (!defined $ref_seq{$chr});

  my $chr_len = $ref_len{$chr};
  my $chr_seq = $ref_seq{$chr};

  if ($cac_strand eq "+") {
    my $up_s = _max(1, $cpos_abs - $slop);
    my $up_e = $cpos_abs - 1;
    my $dn_s = $cpos_abs;
    my $dn_e = _min($chr_len, $cpos_abs + $slop - 1);

    my $up = "";
    $up = substr($chr_seq, $up_s - 1, $up_e - $up_s + 1) if ($up_s <= $up_e);
    my $dn = substr($chr_seq, $dn_s - 1, $dn_e - $dn_s + 1);

    return lc($up) . uc($dn);
  }

  if ($cac_strand eq "-") {
    my $dn_s = _max(1, $cpos_abs - ($slop - 1));
    my $dn_e = $cpos_abs;
    my $up_s = $cpos_abs + 1;
    my $up_e = _min($chr_len, $cpos_abs + $slop);

    my $dn_gen = substr($chr_seq, $dn_s - 1, $dn_e - $dn_s + 1);
    my $dn = revcomp($dn_gen);

    my $up = "";
    if ($up_s <= $up_e) {
      my $up_gen = substr($chr_seq, $up_s - 1, $up_e - $up_s + 1);
      $up = revcomp($up_gen);
    }

    return lc($up) . uc($dn);
  }

  return undef;
}

sub get_upstream_flank_on_cac_strand {
  my ($chr, $cpos_abs, $cac_strand, $L) = @_;
  return undef if (!defined $chr || !defined $cpos_abs || $cpos_abs !~ /^\d+$/);
  return undef if (!defined $ref_seq{$chr});

  my $chr_len = $ref_len{$chr};
  my $chr_seq = $ref_seq{$chr};

  if ($cac_strand eq "+") {
    my $s = _max(1, $cpos_abs - $L);
    my $e = $cpos_abs - 1;
    return "" if ($s > $e);
    return substr($chr_seq, $s - 1, $e - $s + 1);
  } else {
    my $s = $cpos_abs + 1;
    my $e = _min($chr_len, $cpos_abs + $L);
    return "" if ($s > $e);
    my $gen = substr($chr_seq, $s - 1, $e - $s + 1);
    return revcomp($gen);
  }
}

sub find_cac_near_junction {
  my ($chr, $junction, $mode) = @_;

  return (undef, undef, undef) if (!defined $chr || !defined $junction || !defined $mode);
  return (undef, undef, undef) if (!defined $ref_seq{$chr});
  return (undef, undef, undef) if ($junction !~ /^\d+$/);

  my $key = $chr . ":" . $junction . ":" . $mode;
  if (exists $cac_cache{$key}) {
    my $v = $cac_cache{$key};
    return ($v->[0], $v->[1], $v->[2]);
  }

  my $chr_len = $ref_len{$chr};
  my ($region_s, $region_e, $ideal_motif, $cryptic_motif, $choose_rightmost);

  if ($mode eq "UP") {
    $region_s = _max(1, $junction - $search_win);
    $region_e = _min($chr_len, ($junction - 1) + $maxIntrudeIntoAlnBp);

    $ideal_motif = $ideal_heptamer_minus_on_plus;
    $cryptic_motif = $cryptic_trimer_minus_on_plus;
    $choose_rightmost = 1;
  } elsif ($mode eq "DN") {
    $region_s = _max(1, ($junction + 1) - $maxIntrudeIntoAlnBp);
    $region_e = _min($chr_len, $junction + $search_win);

    $ideal_motif = $ideal_heptamer_plus_on_plus;
    $cryptic_motif = $cryptic_trimer_plus_on_plus;
    $choose_rightmost = 0;
  } else {
    return (undef, undef, undef);
  }

  if ($region_s > $region_e) {
    $cac_cache{$key} = [undef, undef, undef];
    return (undef, undef, undef);
  }

  my $subseq = substr($ref_seq{$chr}, $region_s - 1, $region_e - $region_s + 1);

  my $ideal_idx   = $choose_rightmost ? rindex($subseq, $ideal_motif)   : index($subseq, $ideal_motif);
  my $cryptic_idx = $choose_rightmost ? rindex($subseq, $cryptic_motif) : index($subseq, $cryptic_motif);

  my ($ideal_del, $ideal_cpos_abs, $ideal_ok) = (0, 0, 0);
  my ($cryptic_del, $cryptic_cpos_abs, $cryptic_ok) = (0, 0, 0);

  if ($ideal_idx >= 0) {
    if ($mode eq "UP") {
      my $match_end = $region_s + $ideal_idx + length($ideal_motif) - 1;
      $ideal_del = ($junction - 1) - $match_end;   # can be negative
      $ideal_cpos_abs = $match_end;
    } else {
      my $match_start = $region_s + $ideal_idx;
      $ideal_del = $match_start - ($junction + 1); # can be negative
      $ideal_cpos_abs = $match_start;
    }
    $ideal_ok = 1;
  }

  if ($cryptic_idx >= 0) {
    if ($mode eq "UP") {
      my $match_end = $region_s + $cryptic_idx + length($cryptic_motif) - 1;
      $cryptic_del = ($junction - 1) - $match_end; # can be negative
      $cryptic_cpos_abs = $match_end;
    } else {
      my $match_start = $region_s + $cryptic_idx;
      $cryptic_del = $match_start - ($junction + 1); # can be negative
      $cryptic_cpos_abs = $match_start;
    }
    $cryptic_ok = 1;
  }

  if (!$ideal_ok && !$cryptic_ok) {
    $cac_cache{$key} = [undef, undef, undef];
    return (undef, undef, undef);
  }

  my ($cac_strand, $cpos_abs, $raw_del);

  if ($ideal_ok && $cryptic_ok) {
    if (($ideal_del - $cryptic_del) > $ideal_vs_cryptic_dist_thres) {
      $raw_del = $cryptic_del;
      $cpos_abs = $cryptic_cpos_abs;
    } else {
      $raw_del = $ideal_del;
      $cpos_abs = $ideal_cpos_abs;
    }
  } elsif ($ideal_ok) {
    $raw_del = $ideal_del;
    $cpos_abs = $ideal_cpos_abs;
  } else {
    $raw_del = $cryptic_del;
    $cpos_abs = $cryptic_cpos_abs;
  }

  $cac_strand = ($mode eq "UP") ? "-" : "+";

  $cac_cache{$key} = [$cac_strand, $cpos_abs, $raw_del];
  return ($cac_strand, $cpos_abs, $raw_del);
}

sub adjust_junction_and_intrudeRSSlen {
  my ($mode, $cpos_abs, $junction_raw, $raw_del) = @_;

  my $intrudeRSSlen = 0;
  my $junction_adj  = $junction_raw;
  my $del_adj       = $raw_del;

  if (defined $raw_del && $raw_del =~ /^-?\d+$/ && $raw_del < 0) {
    $intrudeRSSlen = -$raw_del;
    $del_adj = 0;
    if ($mode eq "DN") {
      $junction_adj = $cpos_abs - 1;
    } else {
      $junction_adj = $cpos_abs + 1;
    }
  }

  return ($junction_adj, $del_adj, $intrudeRSSlen);
}

sub calc_pins_len_cap {
  my ($seq, $end_side, $B_Qend_adj, $Qstart_adj, $up_flank, $maxPinsLen) = @_;

  return 0 if (!defined $seq || $seq eq "");
  return 0 if (!defined $up_flank);

  $seq = uc($seq);
  $up_flank = uc($up_flank);

  my $Lmax = $maxPinsLen;
  $Lmax = length($up_flank) if (length($up_flank) < $Lmax);
  return 0 if ($Lmax <= 0);

  my $ins_dir = "";

  if ($end_side eq "bait") {
    return 0 if (!defined $B_Qend_adj || $B_Qend_adj !~ /^\d+$/);
    my $start1 = $B_Qend_adj + 1;
    return 0 if ($start1 < 1 || $start1 > length($seq));
    $ins_dir = substr($seq, $start1 - 1, $Lmax);
  } else {
    return 0 if (!defined $Qstart_adj || $Qstart_adj !~ /^\d+$/);
    my $end1 = $Qstart_adj - 1;
    return 0 if ($end1 < 1);
    my $start1 = $end1 - $Lmax + 1;
    return 0 if ($start1 < 1);
    my $frag = substr($seq, $start1 - 1, $Lmax);
    $ins_dir = reverse($frag);
  }

  for (my $L = $Lmax; $L >= 1; $L--) {
    my $insL  = substr($ins_dir, 0, $L);
    my $upSuf = substr($up_flank, length($up_flank) - $L, $L);
    my $pal   = revcomp($upSuf);
    return $L if ($insL eq $pal);
  }

  return 0;
}

sub shrink_pins_to_fit_mid {
  my ($midLen_adj, $P_PinsLen, $B_PinsLen) = @_;

  my $p = defined $P_PinsLen ? int($P_PinsLen) : 0;
  my $b = defined $B_PinsLen ? int($B_PinsLen) : 0;

  return ($p, $b) if (!defined $midLen_adj || $midLen_adj !~ /^\d+$/ || $midLen_adj <= 0);
  return ($p, $b) if (($p + $b) <= $midLen_adj);

  while (($p + $b) > $midLen_adj) {
    if ($p > $b) {
      $p-- if ($p > 0);
      next;
    }
    if ($b > $p) {
      $b-- if ($b > 0);
      next;
    }

    # equal
    if (($p + $b) - $midLen_adj >= 2) {
      $p-- if ($p > 0);
      $b-- if ($b > 0);
    } else {
      # need shrink by 1
      if ($p > 0) { $p--; }
      elsif ($b > 0) { $b--; }
      else { last; }
    }
  }

  $p = 0 if ($p < 0);
  $b = 0 if ($b < 0);
  return ($p, $b);
}

# -------------------------------
# Read TLX headline; determine indices
# -------------------------------
open(my $fh0, "<", $in_tlx) or die "Error: cannot open tlx file $in_tlx for input: $!\n";
my $headline = <$fh0>;
defined $headline or die "Error: empty tlx file $in_tlx\n";
$headline =~ s/[\r\n]+$//;

my @hdr = split(/\t/, $headline, -1);
my %idx;
for (my $i = 0; $i < @hdr; $i++) { $idx{$hdr[$i]} = $i; }

my @need_cols = qw(
  Qname Seq
  B_Qend B_Rname B_Rstart B_Rend B_Strand
  Qstart Rname Rstart Rend Strand
);
for my $c (@need_cols) {
  die "Error: missing required column $c in TLX headline.\n" if (!exists $idx{$c});
}

# -------------------------------
# Infer majority bait CAC coord (optional optimization)
# -------------------------------
my $bait_guess_n = 2000;
if (defined $bait_hint && $bait_hint ne "" && $bait_hint =~ /^\d+$/) {
  $bait_guess_n = int($bait_hint);
  $bait_guess_n = 1 if ($bait_guess_n < 1);
}

my ($maj_b_chr, $maj_b_rs, $maj_b_re, $maj_b_strand, $maj_b_junc) = (undef, undef, undef, 0, undef);
my ($maj_b_cac_coord, $maj_b_cac_strand, $maj_b_cac_abs) = (undef, undef, undef);

if (defined $bait_hint && $bait_hint ne "" && $bait_hint !~ /^\d+$/) {
  my ($cchr, $csignpos) = parse_coord($bait_hint);
  if (defined $cchr && defined $csignpos && $csignpos != 0) {
    $maj_b_cac_strand = ($csignpos < 0) ? "-" : "+";
    $maj_b_cac_abs    = ($csignpos < 0) ? -$csignpos : $csignpos;
    $maj_b_cac_coord  = $bait_hint;
    print STDERR "Note: bait_CAC_coord provided by user: $maj_b_cac_coord\n";
  } else {
    print STDERR "Warning: invalid bait_hint; will try to infer majority bait CAC from TLX.\n";
  }
}

if (!defined $maj_b_cac_coord) {
  my %bait_mode_cnt;
  my $seen = 0;

  while (my $line = <$fh0>) {
    last if ($seen >= $bait_guess_n);
    $line =~ s/[\r\n]+$//;
    next if ($line =~ /^\s*$/);
    my @f = split(/\t/, $line, -1);
    $seen++;

    my $b_chr    = $f[$idx{"B_Rname"}];
    my $b_rs     = $f[$idx{"B_Rstart"}];
    my $b_re     = $f[$idx{"B_Rend"}];
    my $b_strand = strand_sign($f[$idx{"B_Strand"}]);

    next if (!defined $b_chr || $b_chr eq "" || $b_rs !~ /^\d+$/ || $b_re !~ /^\d+$/ || $b_strand == 0);
    my $k = join("\t", $b_chr, $b_rs, $b_re, $b_strand);
    $bait_mode_cnt{$k}++;
  }

  my ($maj_k, $maj_cnt) = ("", 0);
  for my $k (keys %bait_mode_cnt) {
    if ($bait_mode_cnt{$k} > $maj_cnt) { $maj_cnt = $bait_mode_cnt{$k}; $maj_k = $k; }
  }

  if ($maj_k ne "") {
    ($maj_b_chr, $maj_b_rs, $maj_b_re, $maj_b_strand) = split(/\t/, $maj_k);
    $maj_b_junc = ($maj_b_strand > 0) ? $maj_b_re : $maj_b_rs;
    print STDERR "Note: bait majority (first $bait_guess_n): chr=$maj_b_chr B_Rstart=$maj_b_rs B_Rend=$maj_b_re B_Strand=$maj_b_strand (count=$maj_cnt)\n";

    if (defined $maj_b_chr && defined $ref_seq{$maj_b_chr} && defined $maj_b_junc) {
      my $mode = ($maj_b_strand > 0) ? "DN" : "UP";
      my ($cac_strand, $cpos_abs, $raw_del) = find_cac_near_junction($maj_b_chr, $maj_b_junc, $mode);
      if (defined $cpos_abs) {
        $maj_b_cac_strand = $cac_strand;
        $maj_b_cac_abs    = $cpos_abs;
        $maj_b_cac_coord  = format_coord($maj_b_chr, $cac_strand, $cpos_abs);
        print STDERR "Note: inferred majority bait CAC: $maj_b_cac_coord (raw_B_delLen=$raw_del)\n";
      } else {
        print STDERR "Warning: cannot infer majority bait CAC near $maj_b_chr:$maj_b_junc\n";
      }
    }
  } else {
    print STDERR "Warning: cannot determine bait majority from first $bait_guess_n lines.\n";
  }
}

close($fh0);

# -------------------------------
# Second pass: annotate all reads
# -------------------------------
open(my $fh2, "<", $in_tlx) or die "Error: cannot open tlx file $in_tlx for input: $!\n";
my $headline2 = <$fh2>;
$headline2 =~ s/[\r\n]+$//;

print $headline2,
  "\tmidLen",
  "\tB_CAC_coord\tB_CAC_ref_seq\tB_delLen\tB_delLen_adj",
  "\tP_CAC_coord\tP_CAC_ref_seq\tP_delLen\tP_delLen_adj",
  "\tB_Qend_adj\tB_Junction_adj\tQstart_adj\tJunction_adj",
  "\tmidLen_adj\tinsLen",
  "\tP_PelementLen\tB_PelementLen\tnonPinsLen\tMHlen",
  "\n";

my %warn_missing_ref_chr;
my $NR = 1;

while (my $line = <$fh2>) {
  $NR++;
  $line =~ s/[\r\n]+$//;
  next if ($line =~ /^\s*$/);

  my @f = split(/\t/, $line, -1);

  my $Qname    = $f[$idx{"Qname"}];
  my $seq_read = $f[$idx{"Seq"}];

  if (defined $debug_Qname && defined $Qname && $Qname eq $debug_Qname) {
    print STDERR $line . "\n";
  }

  my $B_Qend_raw = $f[$idx{"B_Qend"}];
  my $B_Rname    = $f[$idx{"B_Rname"}];
  my $B_Rstart   = $f[$idx{"B_Rstart"}];
  my $B_Rend     = $f[$idx{"B_Rend"}];
  my $B_Strand   = strand_sign($f[$idx{"B_Strand"}]);

  my $Qstart_raw = $f[$idx{"Qstart"}];
  my $Rname      = $f[$idx{"Rname"}];
  my $Rstart     = $f[$idx{"Rstart"}];
  my $Rend       = $f[$idx{"Rend"}];
  my $Strand     = strand_sign($f[$idx{"Strand"}]);

  if (!exists($ref_seq{$Rname})) {
    if (!exists($warn_missing_ref_chr{$Rname})) {
      print STDERR "Warning: (Line $NR $Qname) Rname $Rname not found in ref.fa $ref_fa\n";
      $warn_missing_ref_chr{$Rname} = 1;
    }
    next;   # skip the records with Rname not in ref.fa
  }
  if (!exists($ref_seq{$B_Rname})) {
    if (!exists($warn_missing_ref_chr{$B_Rname})) {
      print STDERR "Warning: (Line $NR $Qname) B_Rname $B_Rname not found in ref.fa $ref_fa\n";
      $warn_missing_ref_chr{$B_Rname} = 1;
    }
    next;   # skip the records with B_Rname not in ref.fa
  }

  # midLen (raw)
  my $midLen = undef;
  if (defined $Qstart_raw && $Qstart_raw =~ /^\d+$/ && defined $B_Qend_raw && $B_Qend_raw =~ /^\d+$/) {
    $midLen = $Qstart_raw - $B_Qend_raw - 1;
  }

  # ---------- PREY ----------
  my ($P_CAC_coord, $P_CAC_ref_seq) = (undef, undef);
  my ($P_delLen, $P_delLen_adj) = (undef, undef);
  my $Junction_adj = undef;
  my ($P_mode, $P_cac_strand, $P_cpos_abs) = (undef, undef, undef);
  my $P_raw_del = undef;
  my $P_intrudeRSSlen = 0;

  if ($Strand != 0 && defined $ref_seq{$Rname} && defined $Rstart && $Rstart =~ /^\d+$/ && defined $Rend && $Rend =~ /^\d+$/) {
    my $Junction_raw = ($Strand > 0) ? $Rstart : $Rend;
    $P_mode = ($Strand > 0) ? "UP" : "DN";

    ($P_cac_strand, $P_cpos_abs, $P_raw_del) = find_cac_near_junction($Rname, $Junction_raw, $P_mode);

    if (defined $P_cpos_abs) {
      $P_CAC_coord   = format_coord($Rname, $P_cac_strand, $P_cpos_abs);
      $P_CAC_ref_seq = build_cac_seq($Rname, $P_cpos_abs, $P_cac_strand, $CAC_slop_bp);

      $P_delLen = $P_raw_del;

      my ($j_adj, $del_adj, $intrude_len) = adjust_junction_and_intrudeRSSlen($P_mode, $P_cpos_abs, $Junction_raw, $P_raw_del);
      $Junction_adj    = $j_adj;
      $P_delLen_adj    = $del_adj;
      $P_intrudeRSSlen = $intrude_len;
    } else {
      print STDERR "Warning: (Line $NR $Qname) cannot locate prey CAC near $Rname:$Junction_raw (mode=$P_mode)\n";
    }
  }

  # ---------- BAIT ----------
  my ($B_CAC_coord, $B_CAC_ref_seq) = (undef, undef);
  my ($B_delLen, $B_delLen_adj) = (undef, undef);
  my $B_Junction_adj = undef;
  my ($B_mode, $B_cac_strand, $B_cpos_abs) = (undef, undef, undef);
  my $B_raw_del = undef;
  my $B_intrudeRSSlen = 0;

  if ($B_Strand != 0 && defined $ref_seq{$B_Rname} && defined $B_Rstart && $B_Rstart =~ /^\d+$/ && defined $B_Rend && $B_Rend =~ /^\d+$/) {
    my $B_Junction_raw = ($B_Strand > 0) ? $B_Rend : $B_Rstart;
    $B_mode = ($B_Strand > 0) ? "DN" : "UP";

    my $reuse_majority = 0;
    if (defined $maj_b_cac_coord
        && defined $maj_b_chr && defined $maj_b_junc
        && $B_Rname eq $maj_b_chr && $B_Strand == $maj_b_strand
        && $maj_b_junc =~ /^\d+$/ && $B_Junction_raw =~ /^\d+$/) {
      my $d = $B_Junction_raw - $maj_b_junc;
      $d = -$d if ($d < 0);
      $reuse_majority = 1 if ($d <= $maj_reuse_max_dist);
    }

    if ($reuse_majority) {
      if (defined $maj_b_cac_abs && defined $maj_b_cac_strand) {
        $B_cac_strand  = $maj_b_cac_strand;
        $B_cpos_abs    = $maj_b_cac_abs;
        $B_CAC_coord   = $maj_b_cac_coord;
        $B_CAC_ref_seq = build_cac_seq($B_Rname, $B_cpos_abs, $B_cac_strand, $CAC_slop_bp);

        if ($B_mode eq "DN") {
          $B_raw_del = $B_cpos_abs - ($B_Junction_raw + 1);
        } else {
          $B_raw_del = ($B_Junction_raw - 1) - $B_cpos_abs;
        }
      } else {
        print STDERR "Warning: (Line $NR $Qname) majority bait CAC unavailable; cannot reuse for $B_Rname:$B_Junction_raw\n";
        $reuse_majority = 0;
      }
    }

    if (!$reuse_majority) {
      ($B_cac_strand, $B_cpos_abs, $B_raw_del) = find_cac_near_junction($B_Rname, $B_Junction_raw, $B_mode);
      if (defined $B_cpos_abs) {
        $B_CAC_coord   = format_coord($B_Rname, $B_cac_strand, $B_cpos_abs);
        $B_CAC_ref_seq = build_cac_seq($B_Rname, $B_cpos_abs, $B_cac_strand, $CAC_slop_bp);
      }
    }

    if (defined $B_cpos_abs) {
      $B_delLen = $B_raw_del;

      my ($j_adj, $del_adj, $intrude_len) = adjust_junction_and_intrudeRSSlen($B_mode, $B_cpos_abs, $B_Junction_raw, $B_raw_del);
      $B_Junction_adj  = $j_adj;
      $B_delLen_adj    = $del_adj;
      $B_intrudeRSSlen = $intrude_len;
    } else {
      print STDERR "Warning: (Line $NR $Qname) cannot locate bait CAC near $B_Rname:$B_Junction_raw (mode=$B_mode)\n";
    }
  }

  # ---------- Adjust read boundary coords ----------
  my ($B_Qend_adj, $Qstart_adj) = (undef, undef);

  if (defined $B_Qend_raw && $B_Qend_raw =~ /^\d+$/) {
    $B_Qend_adj = $B_Qend_raw - $B_intrudeRSSlen;
    $B_Qend_adj = 0 if ($B_Qend_adj < 0);
  }

  if (defined $Qstart_raw && $Qstart_raw =~ /^\d+$/) {
    $Qstart_adj = $Qstart_raw + $P_intrudeRSSlen;
  }

  # ---------- midLen_adj and insLen ----------
  my ($midLen_adj, $insLen) = (undef, undef);
  if (defined $Qstart_adj && $Qstart_adj =~ /^\d+$/ && defined $B_Qend_adj && $B_Qend_adj =~ /^\d+$/) {
    $midLen_adj = $Qstart_adj - $B_Qend_adj - 1;
    $insLen = ($midLen_adj > 0) ? $midLen_adj : 0;
  }

  # ---------- PinsLen ----------
  my ($P_PinsLen, $B_PinsLen) = (0, 0);

  if (defined $midLen_adj && $midLen_adj =~ /^\d+$/ && $midLen_adj > 0) {
    if (defined $P_delLen_adj && $P_delLen_adj =~ /^\d+$/ && $P_delLen_adj == 0 && defined $P_cpos_abs) {
      my $up = get_upstream_flank_on_cac_strand($Rname, $P_cpos_abs, $P_cac_strand, $maxPinsLen);
      $P_PinsLen = calc_pins_len_cap($seq_read, "prey", $B_Qend_adj, $Qstart_adj, $up, $maxPinsLen);
    }
    if (defined $B_delLen_adj && $B_delLen_adj =~ /^\d+$/ && $B_delLen_adj == 0 && defined $B_cpos_abs) {
      my $up = get_upstream_flank_on_cac_strand($B_Rname, $B_cpos_abs, $B_cac_strand, $maxPinsLen);
      $B_PinsLen = calc_pins_len_cap($seq_read, "bait", $B_Qend_adj, $Qstart_adj, $up, $maxPinsLen);
    }
  }

  # ---------- NinsLen and MHlen ----------
  my $NinsLen = undef;
  my $MHlen = undef;
  if (defined $midLen_adj && $midLen_adj =~ /^-?\d+$/) {
    my $n = $midLen_adj - $P_PinsLen - $B_PinsLen;
    if ($n >= 0){
      $NinsLen = $n;
      $MHlen = 0;
    }else{
      $MHlen = -$n;
      $NinsLen = 0;
    }
  }

  print $line,
    "\t", outv($midLen),

    "\t", outv($B_CAC_coord), "\t", outv($B_CAC_ref_seq), "\t", outv($B_delLen), "\t", outv($B_delLen_adj),
    "\t", outv($P_CAC_coord), "\t", outv($P_CAC_ref_seq), "\t", outv($P_delLen), "\t", outv($P_delLen_adj),

    "\t", outv($B_Qend_adj), "\t", outv($B_Junction_adj), "\t", outv($Qstart_adj), "\t", outv($Junction_adj),

    "\t", outv($midLen_adj), "\t", outv($insLen),

    "\t", outv($P_PinsLen), "\t", outv($B_PinsLen),
    "\t", outv($NinsLen), "\t", outv($MHlen),
    "\n";
}

close($fh2);
