while( <> ) {
chomp;
push @r, [ split "," ];
}
for $i ( 1..scalar( @{ $r[0] } )-
1){#For
each column...
( $min, $ma
x)=(
100000, 0 );      # ... find min and max
for $j ( 0..scalar @r-1 ) {          # ... over all rows.
$min = $r[$j][$i] < $min ? $r[$j][$i] : $min;
$max = $r[$j][$i] > $max ? $r[$j][$i] : $max;
}
for $j ( 0..scalar @r-1 ) {     # Rescale this column in all rows
$r[$j][$i] = ($r[$j][$i] - $min)/($max-$min);
}
}
for $r ( @r ) {
# unless( $r->[3] > 0.6
7 ) { next
; } # Optional filter logic
for $i ( 1..scalar( @{ $r[0] } )-1 ) {
print "$i\t", $r->[$i], "\n";
}
print "\n\n";
}
