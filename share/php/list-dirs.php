<?php
$dirs = array_filter(glob('*'), 'is_dir');
echo json_encode($dirs);
?>
