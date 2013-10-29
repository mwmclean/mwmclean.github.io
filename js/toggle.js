function toggle(showHideDiv, t1, t2, t3, t4, t5, t6) {
	var ele = document.getElementById(showHideDiv);
    if(ele.style.display == "none") {
		ele.style.display = "block";
	}
    var t1 = document.getElementById(t1);
    var t2 = document.getElementById(t2);
    var t3 = document.getElementById(t3);
    var t4 = document.getElementById(t4);
    var t5 = document.getElementById(t5);
    var t6 = document.getElementById(t6);
    t1.style.display = "none";
    t2.style.display = "none";
    t3.style.display = "none";
    t4.style.display = "none";
    t5.style.display = "none";
    t6.style.display = "none";
	
	// change background of current navDiv
	$(".navDiv a").each(function(){
		if($.trim(this.innerHTML.toLowerCase()) == $.trim(showHideDiv)){
			$(this).parent().css( "background-color", "#ffff80" );
			$(this).css( "border-bottom", "#ffff80" );
		}else{
			$(this).parent().css("background-color", "#FFFFFF");
			$(this).css( "border-bottom", "#FFFFFF" );
		}
	});
    var w=window.outerWidth;
    if (showHideDiv == "about") 
    {
        if(w>=1100){
            document.getElementById("thanks").style.display="block";
        }
   //     document.getElementById("thanksTwo").style.display="none";
   //     document.getElementById("showJeff").rel="license";
   //     document.getElementById("showShiny").rel="nofollow";
/*    }
    else if (showHideDiv == "random")
    {
        if(w>=1100){
            document.getElementById("thanksTwo").style.display="block";
        }
        document.getElementById("thanks").style.display="none";
        document.getElementById("showShiny").rel="license";
        document.getElementById("showJeff").rel="nofollow";        */
    }
    else 
    {
        document.getElementById("thanks").style.display="none";
   //     document.getElementById("thanksTwo").style.display="none";
   //     document.getElementById("showJeff").rel="nofollow"; 
   //     document.getElementById("showShiny").rel="nofollow"; 
    }

} 

