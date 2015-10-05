function toggle(folder) {
	//     //display the selected 
	// $(".example-btn").each(function(){
	//     var btn = '#' + divToShow + 'btn';
	//     if(this.id == divToShow){ // || !this.hasClass('collapsed')){
	    var file = '../' + folder + '/index.html';
    // $(btn).trigger("click");
    var divToShow = "ex-" + folder;
    $("#" + divToShow).load(file);
    $("#" + divToShow).toggleClass("in");
    $(".panel-body").each(function(){
    	if(this.id != divToShow && $(this).hasClass("in")){
    	    $(this).toggleClass("in");
    	}
    });
    // $('#' + folder + '-btn').trigger('click');
	// }
	
    // });
    
	
    // change background of btns
    var btn = folder + "-btn";
    $(".example-btn").each(function(){
	if(this.id == btn || $(this).hasClass("shown")){
	    $(this).toggleClass("shown");
	    $(this).blur();  // so bootstrap css for focus'ed btn not used
	}
    });
} 


    	// $(".example-btn").each(function(){
    	//     if(this.id == btn){
	// 	$(this).toggleClass("shown");
    	//       $(this).css( "background-color", "#191919" );
    	//       $(this).css( "color", "white" );
    	//     }else{
    	// 	console.log($(this).css( "background-color"));
    	//     $(this).css( "background-color", "#eceeef" );
    	//     $(this).css( "color", "#191919" );
    	//   }
    	// });
