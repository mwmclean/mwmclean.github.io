function toggleRand(showHideDiv, hide1) {
	var ele = document.getElementById(showHideDiv);
    if(ele.style.display == "none") {
		ele.style.display = "block";
	}
    var hide1 = document.getElementById(hide1);
    hide1.style.display = "none";
	
	if(showHideDiv=='shiny'){
		$('#shinyframe').attr('src', function() {
			return $(this).data('src');
		});
	}
} 

