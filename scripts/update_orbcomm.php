#! /usr/bin/php
<?

$url_path = "http://www.tle.info/data/";
$orbcomm_file = "orbcomm.txt";
$url = $url_path.$orbcomm_file;

echo "Grabbing latest Orbcomm TLE file from $url...\n\n";
echo `wget $url`;

echo "Parsing $orbcomm_file...\n";
$fh = fopen($orbcomm_file, 'r') or die("Can't open file");
$data = fread($fh, filesize($orbcomm_file));
fclose($fh);

$file_cnt = 0;
$tle_header = explode(".",$orbcomm_file);
$tle_header = $tle_header[0];
$new_tle = "$tle_header$file_cnt.tle";

$limit_cnt = 0;
$d_list = split('ORB',$data);
$file_out = fopen($new_tle, 'w');
print "Writing $new_tle...\n";
for ($i=0;$i<sizeof($d_list);$i++) {
	if ($limit_cnt < 21){
		if ($limit_cnt == 0) {
			$out_str = $d_list[$i];
		} else { $out_str = 'ORB'.$d_list[$i];}
		fwrite($file_out, $out_str);
		$limit_cnt++;
	} else {
		fclose($file_out);
		$file_cnt++;
		$limit_cnt = 0;
		$new_tle = "$tle_header$file_cnt.tle";
		print "Writing $new_tle...\n";
		$file_out = fopen($new_tle, 'w');
		$out_str = 'ORB'.$d_list[$i];
		fwrite($file_out, $out_str);
		$limit_cnt++;
	}
}
fclose($file_out);

?>