#!/usr/bin/perl -w
#----------------------------------------------------------------------
use LWP::UserAgent;
#----------------------------------------------------------------------
$ua = LWP::UserAgent->new();
#----------------------------------------------------------------------

if ($#ARGV <= 1)
    {
    print STDERR 'Usage: UserUtil <Host>[:port] <login:password> <[n]{on|off|pulse|status}> ...'."\n";
    exit -1;
    }
($epc, $auth)=splice(@ARGV,0,2);
$base='http://'.$auth.'@'.$epc.'/';

foreach (@ARGV)
{
    $_=lc;
    s/(^[^1-8])/a$1/;
    if (/^([1-8a])on$/)
	{
	RelLink('outlet?'.$1.'=ON');
	}
    elsif (/^([1-8a])off$/)
	{
	RelLink('outlet?'.$1.'=OFF');
	}
    elsif (/^([1-8a])pulse$/)
	{
	RelLink('outlet?'.$1.'=CCL');
	}
    elsif (/^([1-8a])status$/)
	{
	$n=$1;
	defined($response) && ($response->content =~/<a href=outleto/) || RelLink('');
	$content=$response->content;
	while ($content =~ /<a href=outlet\?([1-8])=(on|off)>/ig)
	    {
	    if (($1 eq $n) || ($n eq 'a'))
		{
		if ($2 eq "ON")
		    {print $1," OFF\n";}
		else
		    {print $1," ON\n";}
		}
	    }
	}
    else
	{
	die "Unknown command $_\n";
	}
}

sub RelLink
{
local ($_) = @_;
print STDERR $base.$_,"\n";
$response = $ua->get($base.$_);
$response->is_error() && die $response->status_line;
}

