#!/usr/bin/env perl                                                                                                                         

use strict;
use CGI;

my $cgi = new CGI;
print $cgi->header();
my $jobid = $cgi->param("jobid");

#$pngfileD = "1234972955";
#my $hoge = `ls -l $pngfileD`;
#my @hoge = split(/\s/, $hoge);
#my $hoge = `ls -ld $pngfileD`;
#my $hoge = `ls -ld $pngfileD`;
#my @hoge = split(/\s/, $hoge);
#$jobid = "134019832617892";
#$jobid = "134020250130494";


my $xml = "./Results/".$jobid.".xml";

if(`grep "/Items" $xml`){
    print qq(
             <html>
             <head>
             <meta http-equiv="Refresh" content="1;URL=http://ws.g-language.org/TALEN/index.cgi?jobid=$jobid">
	     );
}else{
    print qq(
             <html>
             <head>
             <meta http-equiv="Refresh" content="6;URL=http://ws.g-language.org/TALEN/running.cgi?jobid=$jobid">
<style type="text/css">
	     body {
	       color: #000;
		 font-family: Arial, Helvetica, sans-serif;
		 font-size: 10pt;
	     }
	     html, body {
	       margin: 0;
	       padding: 0;
	     }
	     /*\*/
	     html, body, .outer, .middle {
	       height: 100%;
	     }
	     .outer {
	       margin: 0 auto;
	       display: table;
	     }
	     .middle {
	       display: table-cell;
		 vertical-align: middle;
	     }
	     /**/
	     .inner {
	       width: 45em;
		 min-width: 700px;
	       padding: 0;
	     }
	     #box {
	   border: 1px solid #333;
	     background-color: #eee;
	 }
    h1 {
        background: url("http://ws.g-language.org/TALEN/images/results_talen_p.jpg") no-repeat top left #EEEEEE;
	height: 70px;
	vertical-align: middle;
	font-size: 10pt;
	font-weight: bold;
	text-align: right;
	padding: 20px 10px 0 0;
	margin: 0;
    }
    h1 span {
	line-height: normal;
	font-weight: normal;
	font-size: 10pt;
      display: block;
    }

    #jobInfo {
  margin: 10px 0;
    text-align: center;
    font-size: 10pt;
    background-color: #FFF;
  padding: 5px 0;
}

#jobInfo h2 {
font-size: 11pt;
padding: 5px 0 10px 0;
margin: 0;
}

h3 {
    text-align: center;
font-size: 10pt;
}

ul li {
    padding-bottom: 10px;
}

li ul li {
    padding-bottom: 3px;
}

#copyright {
text-align: center;
color: #999;
  padding: 5px;
}
</style>
    <!--[if IE]>
    <style>
    .middle {
	text-align: center;
    }
.inner, .ieSpan {
    vertical-align: middle;
}
.inner {
  display: inline;
  _height: 0;
  zoom: 1;
    text-align: left;
}
.ieSpan {
  height: 100%;
  zoom: 1;
}
</style>
    <![endif]-->
    <!--[if IE 8]>
    <style>
    .inner {
      display: inline-block;
      width: auto;
    }
</style>
             <script language="JavaScript">
             <!--
             function reLoad(){
                 window.location='http://ws.g-language.org/TALEN/running.cgi'
                     setTimeout ("reLoad()",10000);
             }
             //-->
             </script>
             </head>
             <body onload="setTimeout('reLoad()',10000)">
             <div class="outer">
	     <div class="middle">
	     <div class="inner">
	     <div id="box">
	     <h1>Your job [$jobid] is currently running... <span>please be patient</span></h1>
	     <div id="jobInfo">
             <h2>The TALEN and primer data will appear in this window.</h2>
	     </div>
             <center><img src="./images/loading.gif" width="16" height="16"/></center>
	     <h3>Please wait.</h3>
	     </div>
	     <div id="copyright">
             &copy; 2012 Institute for Advanced Biosciences, Keio University
	     </div>
	     </div>
	     <span class="ieSpan"></span>
	     </div>
	     </div>
  </body>
</html>
             );
}
