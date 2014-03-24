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
			$(this).parent().css( "background-color", "#ffff80" );
			$(this).css( "border-bottom", "#ffff80" );
		}else{
			$(this).parent().css("background-color", "#FFFFFF");
			$(this).css( "border-bottom", "#FFFFFF" );
		}
		// handle cv div specially
		if ($.trim(this.innerHTML.toLowerCase()) == "html" && divToShow == "cv"){
		   
			$( "div.navDiv.nav-menu" ).css( "background-color", "#ffff80" );
		}else if ($.trim(this.innerHTML.toLowerCase()) == "html"){
			$( "div.navDiv.nav-menu" ).css("background-color", "#FFFFFF");
		}
	});
	
	//display thank you to Jeff if browser window at least as big as specified size
    var w=window.outerWidth;
    if (divToShow == "about") 
    {
        if(w>=1100){
            document.getElementById("thanks").style.display="block";
        }
    }
    else 
    {
        document.getElementById("thanks").style.display="none";
    }

} 

