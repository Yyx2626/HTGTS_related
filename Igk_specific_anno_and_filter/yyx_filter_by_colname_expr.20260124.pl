#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Version: 0.1.0 (2026-01-24)
Authors: Adam Yongxin Ye @ BCH & ChatGPT
';

my $usage = "Usage: perl $0 <input.tlx> <output.pass.tlx> <output.filtered_out.tlx>
    <na_filtered_out(suggest:1)> <filter_pass_expr> [filter_pass_expr_2] ...

Examples of <filter_pass_expr>:
  B_delLen_adj <= 7 && P_delLen_adj <= 7 && insLen <= 2
  (B_delLen_adj + P_delLen_adj) <= 10 and insLen == 0
  B_delLen_adj <= 5 or P_delLen_adj <= 5 && insLen <= 2
  F[5] <= 10 && insLen <= 2
  \$6 <= 10 && \$F[12] <= 2

Notes:
  - If <filter_pass_expr> contains spaces or special characters, it often needs to be quoted
    in your shell (single quotes recommended when using \$N).
  - You may provide multiple filter_pass_expr arguments:
      <expr1> [expr2] [expr3] ...
    They will be combined as:
      expr1 AND expr2 AND expr3 ...
  - Supported bareword operators (lowercase):
      and/or/not   (logical operators)
      eq/ne/lt/le/gt/ge   (string comparisons)
    You may also use Perl symbolic operators such as:
      && || !   (logical operators)
      == != < <= > >=   (numeric comparisons)
      + - * / ( )   (numeric calculations)
  - Field indexing shortcuts:
      F[i] means \$F[i] (0-based index, Perl-style).
      \$N follows awk convention (1-based fields):
        \$1 means the 1st column -> \$F[0]
        \$6 means the 6th column -> \$F[5]
    If you use \$N or \$F[i] in your shell command, use single quotes around <filter_pass_expr>
    (or escape the \$) to prevent the shell from expanding it.
  - Bracket usage:
      In raw <filter_pass_expr>, brackets are only allowed as F[<int>] or \$F[<int>].
      After expansion, brackets are only allowed as \$F[<int>] or \$F[\$f2i{\"<colname>\"}].
  - na_filtered_out:
      When na_filtered_out != 0, if eval triggers warnings like:
        Argument is not numeric
      the record will be written to output.filtered_out.tlx (and the warning still goes to STDERR).
".$version;

sub _die_with_expr {
    my ($msg, $expr_raw, $expr_now) = @_;
    $expr_now = "" if !defined $expr_now;
    die "$msg\n  expr_raw = '$expr_raw'\n  expr_now = '$expr_now'\n";
}

sub assert_only_allowed_brackets {
    my ($expr, $expr_raw, $phase, $f2i_href, $ncol) = @_;

    my $tmp = $expr;
    my $f2i = $f2i_href // {};

    if ($phase eq 'raw') {
        # allow: F[<int>] or $F[<int>]
        $tmp =~ s/\bF\s*\[\s*\d+\s*\]//g;
        $tmp =~ s/\$F\s*\[\s*\d+\s*\]//g;
    }
    elsif ($phase eq 'expanded') {

        # (A) check $F[<int>] index bounds (if $ncol provided)
        if (defined $ncol) {
            while ($tmp =~ /\$F\s*\[\s*(\d+)\s*\]/g) {
                my $idx = $1;
                if ($idx >= $ncol) {
                    _die_with_expr("Error: \$F[$idx] is out of range (ncol=$ncol).", $expr_raw, $expr);
                }
            }
        }

        # (B) check $f2i{\"col\"} keys are valid
        while ($tmp =~ /\$F\s*\[\s*\$f2i\{\"([^\"]+)\"\}\s*\]/g) {
            my $col = $1;
            if (!exists $f2i->{$col}) {
                _die_with_expr("Error: unknown column name used in expression: '$col'.", $expr_raw, $expr);
            }
        }

        # remove allowed: $F[<int>] or $F[$f2i{"..."}]
        $tmp =~ s/\$F\s*\[\s*\d+\s*\]//g;
        $tmp =~ s/\$F\s*\[\s*\$f2i\{\"[^\"]+\"\}\s*\]//g;
    }
    else {
        die "Error: internal error: unknown bracket-check phase '$phase'\n";
    }

    if ($tmp =~ /[\[\]]/) {
        _die_with_expr("Error: illegal bracket usage detected ($phase phase).", $expr_raw, $expr);
    }
    return;
}

sub compile_expr {
    my ($expr_raw, $f2i_href, $ncol) = @_;
    my %f2i = %{$f2i_href};

    my $expr = $expr_raw;

    # Block obviously dangerous chars in eval context (do NOT block brackets here)
    # Keep {} blocked in raw phase; expanded phase will contain $f2i{"col"}.
    if ($expr =~ /[;\`{}\\]/) {
        _die_with_expr("Error: unsafe character(s) detected in raw expression.", $expr_raw, $expr_raw);
    }

    # Bracket usage in raw expression must be restricted to F[<int>] or $F[<int>]
    assert_only_allowed_brackets($expr, $expr_raw, "raw", undef, $ncol);

    # Allow shorthand:
    #   F[5]  -> $F[5]   (Perl-style, 0-based)
    #   $6    -> $F[5]   (awk-style, 1-based)
    $expr =~ s/\bF\s*\[\s*(\d+)\s*\]/\$F\[$1\]/g;

    $expr =~ s/\$(\d+)\b/
        my $awk = $1;
        my $idx = $awk - 1;
        die "Error: \$$awk is invalid (awk fields start from \$1)\n" if $idx < 0;
        if (defined $ncol && $idx >= $ncol) {
            die "Error: \$$awk refers to column $awk but ncol=$ncol\n";
        }
        "\$F[$idx]";
    /eg;

    # Replace column names using placeholder strategy to avoid re-replacement.
    # Use a rare placeholder prefix and salt it with PID to avoid collisions.
    my $salt = $$;
    my @keys = sort { length($b) <=> length($a) } keys %f2i;

    my %ph2key;
    my $k = 0;
    for my $key (@keys) {
        my $ph = "__YyX_cOl_${salt}_$k" . "__";
        $k++;

        # Word-boundary token match; assumes headers are mostly [A-Za-z0-9_]+
        my $n = ($expr =~ s/\b\Q$key\E\b/$ph/g);
        if ($n) {
            $ph2key{$ph} = $key;
        }
    }

    for my $ph (keys %ph2key) {
        my $key = $ph2key{$ph};
        $expr =~ s/\Q$ph\E/\$F[\$f2i{\"$key\"}]/g;
    }

    # Bracket usage after expansion must be restricted + validate f2i keys and $F[] bounds
    assert_only_allowed_brackets($expr, $expr_raw, "expanded", $f2i_href, $ncol);

    # Validate: only allow:
    #   - $F[<int>]
    #   - $F[$f2i{"..."}]
    #   - bareword ops: and/or/not/eq/ne/lt/le/gt/ge
    # Everything else that looks like a bareword token is illegal.
    my $check = $expr;

    $check =~ s/\$F\[\d+\]//g;
    $check =~ s/\$F\[\$f2i\{\"[^\"]+\"\}\]//g;
    $check =~ s/\b(?:and|or|not|eq|ne|lt|le|gt|ge)\b//ig;

    my %bad;
    while ($check =~ /([A-Za-z_][A-Za-z0-9_]*)/g) {
        $bad{$1} = 1;
    }
    if (%bad) {
        my @bad = sort keys %bad;
        _die_with_expr("Error: illegal bareword token(s) in expression code: " . join(", ", @bad), $expr_raw, $expr);
    }

    # Extra guard (expanded): do NOT forbid { } because $f2i{"col"} is expected.
    # Still forbid a few characters we really don't want in eval.
    if ($expr =~ /[;\`\\:]/) {
        _die_with_expr("Error: unsafe character(s) detected in expanded expression.", $expr_raw, $expr);
    }

    return $expr;
}

sub is_numeric_warning {
    my ($w) = @_;
    return 1 if ($w =~ /isn't numeric in numeric/i);
    return 0;
}

# ---------------- main ----------------

if (@ARGV < 5) {
    die $usage;
}

my $input_filename = shift(@ARGV);
my $output_pass_filename = shift(@ARGV);
my $output_filtered_filename = shift(@ARGV);

# Optional na_filtered_out (default:1).
# If next arg is 0 or 1, treat it as na_filtered_out; otherwise keep default=1.
my $na_filtered_out = shift(@ARGV);

my @expr_parts = @ARGV;
if (!@expr_parts) {
    die $usage;
}

# Combine multiple expr parts with AND (Perl: "and")
my $expr_raw = join(" and ", map { "( $_ )" } @expr_parts);

print STDERR "[Filter Expression] na_filtered_out = $na_filtered_out\n";
print STDERR "[Filter Expression] raw = '$expr_raw'\n";

open IN, "<", $input_filename
  or die "Error: cannot open '$input_filename' for input: $!\n";
open OUT_PASS, ">", $output_pass_filename
  or die "Error: cannot open pass '$output_pass_filename' for output: $!\n";
open OUT_FAIL, ">", $output_filtered_filename
  or die "Error: cannot open filtered '$output_filtered_filename' for output: $!\n";

my %f2i;
my $NR = 0;
my $expr_code;
my $header_ncol;

while (my $line = <IN>) {
    $NR++;
    $line =~ s/[\r\n]+$//;

    # Skip empty/blank lines
    next if $line =~ /^\s*$/;

    my @F = split /\t/, $line, -1;

    if ($NR == 1) {
        $header_ncol = scalar(@F);

        for (my $i = 0; $i < @F; $i++) {
            $f2i{$F[$i]} = $i;
        }

        print OUT_PASS join("\t", @F), "\n";
        print OUT_FAIL join("\t", @F), "\n";

        $expr_code = compile_expr($expr_raw, \%f2i, $header_ncol);
        print STDERR "[Filter Expression] code = '$expr_code'\n";
        next;
    }

    # If the data line has fewer columns than header, we keep it as-is.
    # It may lead to undef values; na_filtered_out can help catch numeric warnings.
    my $pass = 0;
    my $saw_numeric_warn = 0;

    {
        my @warn_buf;
        local $SIG{__WARN__} = sub {
            my ($w) = @_;
            push @warn_buf, $w;

            # Always print warnings to STDERR (preserve existing behavior)
            print STDERR $w;

            if ($na_filtered_out && is_numeric_warning($w)) {
                $saw_numeric_warn = 1;
            }
        };

        ## no critic
        $pass = eval $expr_code;
        ## use critic

        if ($@) {
            die "Error: (Line $NR) error in evaluating expression code: $@\n"
              . "Expression: $expr_code\n";
        }
    }

    # If numeric warning occurred, force to filtered_out when na_filtered_out=1
    if ($na_filtered_out && $saw_numeric_warn) {
        $pass = 0;
    }

    if ($pass) {
        print OUT_PASS join("\t", @F), "\n";
    } else {
        print OUT_FAIL join("\t", @F), "\n";
    }
}

close IN;
close OUT_PASS;
close OUT_FAIL;

0;
