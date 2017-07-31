#
#	CAM_namelist.pm			Erik Kluzek
#
#	Perl module to create, parse and check the namelists needed
#	for CAM.
#
#------------------------------------------------------------------------
#
#	Description of methods:
#
#	new ----------------------- Constructor
#	set_namelist -------------- Set specifics of namelist (which land-model to use etc)
#	build --------------------- Create the namelist
#	interactive --------------- Interactively change the values in the namelist
#	parse --------------------- Parse a previous namelist.
#	print --------------------- Print contents of namelist to terminal.
#
#	$Id: CAM_namelist.pm,v 1.19.4.12.18.1 2004/05/24 16:07:45 mvr Exp $
#
use strict;
#use diagnostics;
use Cwd;

package CAM_namelist;
use CAM_config;
use camexp;
use clm2exp;

sub new {
#
# Constructor
#
  my $class = shift;
  my $optsref = shift;

  my $self = {};
  my %opts = %$optsref;
  $self->{'ATMNL'} = undef;                       # Atmospheric model namelist object
  $self->{'LNDNL'} = undef;                       # Land-model namelist object
  $self->{'PARSE_FILE'} = $opts{'infile'};        # User supplied input namelist file
  if ( $opts{namelist} ) {                        # namelist settings from the command-line
    my $filename = ".tmpfile.";
    open( TMPNL, ">$filename" ) || die "ERROR:: Can not open temp file : $filename\n";
    my $namelist = $opts{namelist};
    print TMPNL "$namelist\n";
    close( TMPNL );
    $self->{'PARSE_NLLINE'} = $filename;
  } else {
    $self->{'PARSE_NLLINE'} = undef;              # Temporary-Filename of command-line 
  }

  $self->{'INTERACTIVE'} = $opts{interactive};    # Whether to use interactive prompting or not

  $self->{'NAMELIST'} = $opts{out};               # Output filename of namelist

  $self->{'CFG_OBJ'} = CAM_config->new;           # Configuration variables
  $self->{config_cache_file} = $opts{config};     # config_cache.xml filename

  $self->{LANDMODEL} = undef;                     # Land-model to use

  $self->{'optsref'} = $optsref;

  $self->{'printlev'} = $opts{'printlev'};   # Print level

  my $eol = "\n";
  if ($opts{'interactive'}) { $eol = "\n\n"; }
  $self->{eol} = $eol;       # End of line

  bless( $self, $class );
  return( $self );
}

#============================================================================

sub set_namelists {
#
# Create the namelist objects needed
#
  my $self = shift;

  my $class = ref($self);
  my $nm = "$class\:\:set_namelists";
  my $interactive = $self->{INTERACTIVE};
  #
  # Read in the config_cache.xml file
  #
  my $file = $self->{config_cache_file};
  $self->{'CFG_OBJ'}->read_config_cache( $file );
  $self->{'LANDMODEL'} = $self->{'CFG_OBJ'}->cfg("LANDMODEL");

  #
  # Create CAM and land model namelists
  #
  my $NAMELIST = $self->{'NAMELIST'};  # Output filename

  # CAM namelist
  $self->{'ATMNL'} = camexp->new( $self->{'optsref'}, $self->{'CFG_OBJ'} );
  #
  # Land-model namelist (either clm2)
  #
  my $LANDMODEL = $self->{'LANDMODEL'};
  if ($self->{LANDMODEL} eq 'clm2') {
     $self->{'LNDNL'} = clm2exp->new( $self->{'optsref'}, $self->{'CFG_OBJ'} );
  } else {
    die "ERROR($nm): LANDMODEL = " . $self->{LANDMODEL} . ".  Should be clm2\n";
  }    
}

#============================================================================

sub interactive {
#
# Ask for additional settings interactively
#
  my $self = shift;

  $self->{ATMNL}->change;
  $self->{LNDNL}->change;
}

#============================================================================

sub print {
#
# Print out the resulting namelist
#
  my $self = shift;

  $self->{ATMNL}->print;
  $self->{LNDNL}->print;
}

#============================================================================

sub parse {
#
# Parse a namelist
#
  my $self = shift;
  my $file = shift;

  $self->{ATMNL}->parse( $file );
  $self->{LNDNL}->parse( $file );
}

#============================================================================

sub build {
#
# Build the namelist
#
  my $self = shift;

  my $optsref = $self->{'optsref'};

  # Get values from user specified namelist file
  if ( defined( $self->{'PARSE_FILE'} ) ) {
    $self->parse( $self->{'PARSE_FILE'});
  }

  # Get values from namelist specified on the command-line (these values
  # will overwrite the ones read from the namelist file).
  if ( defined( $self->{'PARSE_NLLINE'} ) ) {
    $self->parse( $self->{'PARSE_NLLINE'});
    my $cmd = "/bin/rm -rf " . $self->{'PARSE_NLLINE'};
    `$cmd`;
  }

  # set all keys to lower case (this should be done where the keys are set)
  $self->{'ATMNL'}->convert_case;
  $self->{'LNDNL'}->convert_case;

  # Run type is set in CAM namelist.  Land model namelist only needs to know
  # the run type to check that "nrevsn" is set for a branch simulation.
  my %settings;
  $settings{RUNTYPE} = $$optsref{runtype};
  $settings{ncdata_vers} = $$optsref{ncdata_vers};

  # Set output values according to precedence of various input values.  Ensure
  # that valid output namelists are produced.
  $self->{'ATMNL'}->set_output_values(%settings);
  $self->{'LNDNL'}->set_output_values(%settings);

  # Allow for additional settings interactively
  if ( defined($self->{'INTERACTIVE'}) && $self->{'INTERACTIVE'} ) {
    $self->interactive;
  }

  # Write the final namelist to an output file.
  my $NAMELIST = $self->{NAMELIST};
  my $eol = $self->{eol};
  if ( $self->{'printlev'} ) { print "Write out namelist to: $NAMELIST $eol"; }
  $self->{'ATMNL'}->Write;
  $self->{'LNDNL'}->Write('Append');
}

1   # to make use or require happy
