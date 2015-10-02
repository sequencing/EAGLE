#!/usr/bin/env perl

if (scalar(@ARGV) == 0) {
  print STDERR "Usage:\n";
  print STDERR "cd <doxygen's html/dot directory>\n";
  print STDERR "mkdir dot2\n";
  print STDERR "for i in *dot; do echo \$i; (dot \$i > dot2/\$i) ; done\n";
  print STDERR "doxygen_dot_concat.pl dot2/*__coll__graph.dot > all.dot1\n";
exit(0);
}

my %knownLabels = ();
my %alreadyPlottedNamedEdges = ();
my %extraConnectionsNeeded = ();
my %nodesToDelete = ();


print "strict digraph G\n";
print "{\n";
#print "ratio=compress\n";
#print "size=\"29.7,21\"\n";
#print "  bgcolor=\"transparent\";\n";
print "  edge [fontname=\"FreeSans.ttf\",fontsize=10,labelfontname=\"FreeSans.ttf\",labelfontsize=10];\n";
print "  node [fontname=\"FreeSans.ttf\",fontsize=10,shape=record];\n";


foreach $fileNum (0 .. $#ARGV) {
  my $filename = $ARGV[$fileNum];
  my %duplicatedLocalNames = ();
  my %deletedNodes = ();

  print STDERR "Opening $filename\n";

  open( INP, "$filename" ) or die "Error: Cannot open $inputFile";
  while (my $line = <INP>) {

    while ($line =~ /\\\n/) {
      my $cont = <INP>;
      $line =~ s/\\\n//;
      $line = $line . $cont;
    }
#print "line=$line\n";

    ($line =~ /Node/) or next;

    if ($line =~ /Node([0-9]+) \[label=(.+), height/) {
      print STDERR "This is a node: $1, $2\n";
      my $nodeNum = $1;
      my $label = $2;

      if ($label eq "\"boost::noncopyable\"") {
        # Delete this node and any connected edges
        print STDERR "deleting this node\n";
        $deletedNodes{"Node${nodeNum}"} = 1;
      }
      else {
        if (defined $knownLabels{$label} and $line !~ / color=grey75/) {
          print STDERR "already known\n";
          $duplicatedLocalNames{"Node${nodeNum}"} = $knownLabels{$label};
          #      print STDERR "---\n";
          #      print STDERR "hash=\n";
          #      while (my ($key,$value) = each(%duplicatedLocalNames)) {
          #        print STDERR "$key => $value\n";
          #      }
          #      print STDERR "---\n";
        } else {
          print STDERR "not known yet\n";
          $knownLabels{$label} = "Node${nodeNum}_${fileNum}";
          
          my $replacementNodeName = $knownLabels{$label};
          $line =~ s/Node[0-9]+ /$replacementNodeName /g;
          print $line;

          # Check if we need to manually connect this node to its parent,
          # when doxygen doesn't automatically does it
          # e.g "vector<MyClass>", pair<...>, ...
          if ($label =~ /vector\\< *(eagle[^ ]*) *\\>/) {
            print STDERR "In pass2: Connect $replacementNodeName to parent ($1)\n";
            $extraConnectionsNeeded{$replacementNodeName} = $1;
          }

          if ($label =~ /deque\\< *([^ ]*) *\\>/) {
            print STDERR "In pass2: Connect $replacementNodeName to parent ($1)\n";
            $extraConnectionsNeeded{$replacementNodeName} = $1;
          }

          if ($label =~ /pair\\< *([^ ,]*)[, ].*\\>/) {
            print STDERR "In pass2: Connect $replacementNodeName to parent ($1)\n";
            $extraConnectionsNeeded{$replacementNodeName} = $1;
          }

          if ($label =~ /\"(eagle::[^ \\<]*) *\\<.*::iterator/) {
            my $dest = "$1\\< Iterator \\>";
            print STDERR "In pass2: Connect $replacementNodeName to parent ($dest)\n";
            $extraConnectionsNeeded{$replacementNodeName . " special"} = $dest;
          }
        }
      }
    }
    elsif ($line =~ /Node([0-9]+) \-> Node([0-9]+) \[/) {
      # This is an edge
      my ($nodeNum1, $nodeNum2) = ($1, $2);
      my $nodeName2 = "Node${nodeNum2}_${fileNum}";
      print STDERR "This is an edge: $nodeNum1, $nodeNum2, $label ($line)\n";

      my $doNotPrint = 0;
      while (my ($key,$value) = each(%deletedNodes)) {
        print STDERR "checking if deleting $key \n";
        if ($line =~ /$key /) {
          print STDERR "deleting $key \n";
          $doNotPrint = 1;
        }
      }

      while (my ($key,$value) = each(%duplicatedLocalNames)) {
        print STDERR "replacing $key -> $value \n";
        if ($line =~ s/$key /$value /g) {
          print STDERR "nodeName2=$nodeName2\n";
          $nodeName2 =~ s/^${key}_${fileNum}$/$value/g;
          print STDERR "nodeName2=$nodeName2\n";
        }
      }

      if ($line =~ /, label=(.+), pos=/) {
        my $label = $1;
        print STDERR "label=$label\n";

        while (my ($key,$value) = each(%alreadyPlottedNamedEdges)) {
          print STDERR "checking ${nodeName2}_${label} against already plotted edge $key\n";
          if ("${nodeName2}_${label}" eq "${key}") {
            print STDERR "yes, already plotted edge $key\n";
            $doNotPrint = 1;
            print STDERR "Later delete Node${nodeNum1}_${fileNum} - it should appear disconnected from the rest\n";
            $nodesToDelete{"Node${nodeNum1}_${fileNum}"} = 1;
            break;
          }
        }

        if (! $doNotPrint) {
          print STDERR "adding ${nodeName2}_${label} to list of alreadyPlottedNamedEdges\n";
          $alreadyPlottedNamedEdges{"${nodeName2}_${label}"} = 1;
        }
      }

      $line =~ s/(Node[0-9]+) /\1_$fileNum /g;
      if (! $doNotPrint) {
        print $line;
      }
    }
    else {
      die "Unrecognised line format: $line";
    }
  }
}


  # Pass 2
  print STDERR "PASS 2\n";
  while (my ($key,$value) = each(%extraConnectionsNeeded)) {
    print STDERR "Adding connection $key -> $value \n";
    $key =~ s/ special//;
    $extraConnectionsNeeded{$replacementNodeName} = $1;

    my $destNodeLabel;
    $destNodeLabel = $value;
    if (defined $knownLabels{$destNodeLabel}) {
      my $destNode = $knownLabels{$destNodeLabel};
      print STDERR "OK found: $destNode ($destNodeLabel) -> $key [dir=back]\n";
      print "$destNode -> $key [dir=back]\n";
      next;
    }

    $destNodeLabel = "\"" . $value . "\"";
    if (defined $knownLabels{$destNodeLabel}) {
      my $destNode = $knownLabels{$destNodeLabel};
      print STDERR "OK found: $destNode ($destNodeLabel) -> $key [dir=back]\n";
      print "$destNode -> $key [dir=back]\n";
      next;
    }

    $destNodeLabel = "\"eagle::model::" . $value . "\"";
    if (defined $knownLabels{$destNodeLabel}) {
      my $destNode = $knownLabels{$destNodeLabel};
      print STDERR "OK found: $destNode ($destNodeLabel) -> $key [dir=back]\n";
      print "$destNode -> $key [dir=back]\n";
      next;
    }

    $destNodeLabel = "\"eagle::genome::" . $value . "\"";
    if (defined $knownLabels{$destNodeLabel}) {
      my $destNode = $knownLabels{$destNodeLabel};
      print STDERR "OK found: $destNode ($destNodeLabel) -> $key [dir=back]\n";
      print "$destNode -> $key [dir=back]\n";
      next;
    }

    print STDERR "No match\n";
  }


print "}\n";


#print STDERR "* Nodes to delete:\ngrep -v ";
print STDERR "\n\n";
print STDERR "* Now you can run:\n";
print STDERR "grep -v ";
while (my ($key,$value) = each(%nodesToDelete)) {
  print STDERR "-e $key ";
}
print STDERR "all.dot1 > all.dot2\n\n";

print STDERR "* Then:\n";
print STDERR "dot all.dot2 | ccomps -x | gvpack > all.dot3\n\n";
print STDERR "* And:\n";
print STDERR "neato -n2 -Tpng all.dot3 > all.png\n";
print STDERR "neato -n2 -Tpdf all.dot3 > all.pdf\n";
