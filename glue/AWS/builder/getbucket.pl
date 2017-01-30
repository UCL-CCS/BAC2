#! /usr/bin/perl

use LWP::Simple;

#get the bucket id
$userdata = get("http://169.254.169.254/latest/user-data");

#print $userdata."\n";

$userdata =~ /BUCKETID: (.*)/;

$bucketid=$1;

print $bucketid;
