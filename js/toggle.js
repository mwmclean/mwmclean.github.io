function toggle(divToShow) {
    //display the selected 
    $(".contentDiv").each(function(){
	if(this.id == divToShow){
		$(this).css( "display", "block" );
	}else{
		$(this).css( "display", "none" );
	}
    });	

    // change background of current navDiv
    $(".navDiv a").each(function(){
	if($.trim(this.innerHTML.toLowerCase()) == divToShow){
	    $(this).parent().css( "background-color", "#FFFFFF" );
	    $(this).css( "text-decoration", "underline" );
	}else{
	    $(this).parent().css("background-color", "#FFFFFF");
	    $(this).css( "text-decoration", "none" );
	}
	// handle cv div specially
	if ($.trim(this.innerHTML.toLowerCase()) == "html" && divToShow == "cv"){
	    $( "div.navDiv.nav-menu" ).css( "text-decoration", "underline" );
	}else if ($.trim(this.innerHTML.toLowerCase()) == "html"){
	    $( "div.navDiv.nav-menu" ).css("text-decoration", "none");
	}
    });	
} 

