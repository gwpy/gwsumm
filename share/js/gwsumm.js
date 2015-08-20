/*
 * Copyright (C) Duncan Macleod (2013)
 *
 * This file is part of GWSumm
 *
 * GWSumm is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GWSumm is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GWSumm.  If not, see <http://www.gnu.org/licenses/>
 */

/* ------------------------------------------------------------------------- */
/* Calendar links                                                            */
var re_dayurl = new RegExp('day\/\\d{8}\/')
var re_monthurl = new RegExp('month\/\\d{6}\/')
var re_yearurl = new RegExp('year\/\\d{4}\/')

function findDateFormat() {
    var url = window.location.href;
    if (re_dayurl.test(url)) {
        return 'day';
    }
    else if (re_monthurl.test(url)) {
        return 'month';
    }
    else if (re_yearurl.test(url)) {
        return 'year';
    }
    return undefined;
}

function getPageDate() {
    var dateformat = findDateFormat();
    if (dateformat == 'day') {
        var datestring = String(re_dayurl.exec(
                                    window.location.href)).split('/')[1];
        return moment(datestring, 'YYYYMMDD');
    }
    else if (dateformat == 'month') {
        var datestring = String(re_monthurl.exec(
                                    window.location.href)).split('/')[1];
        return moment(datestring, 'YYYYMM');
    }
    else if (dateformat == 'year') {
        var datestring = String(re_yearurl.exec(
                                    window.location.href)).split('/')[1];
        return moment(datestring, 'YYYY');
    }
    return undefined;
}

function stepDate(step) {
    var dateformat = findDateFormat();
    if (!dateformat) { return; }
    var date = getPageDate();
    var newdate = date.add(dateformat, step);
    if (dateformat == 'day') {
        move_to_date({type: 'moveDate', date: newdate});
    }
    else if (dateformat == 'month') {
        move_to_date({type: 'moveMonth', date: newdate});
    }
    else if (dateformat == 'year') {
        move_to_date({type: 'moveYear', date: newdate});
    }
}

// Load a given 'state'
$.fn.load_state = function loadState(page) {
    if ($(this).attr('id') == undefined) {
        return
    }
    $('#content').load(page);
    $(this).set_state();
}

// Set a given state in the menu
$.fn.set_state = function setState() {
    if ($(this).attr('id') == undefined) {
        return;
    }
    var id = $(this).attr('id').substring(6);
    $('#states').html($(this).attr('title') + ' <b class="caret"></b>');
    $('.state').removeClass('open');
    $(this).toggleClass('open');
    window.location.hash = '#' + id;
    $('a.ifo-switch').each(function() {
        var oldurl = $(this).attr('href');
        var oldhash = oldurl.split('#')[1];
        $(this).attr('href', oldurl.replace(oldhash, id));
    });
}

// Move to the date selected
function move_to_date(ev) {
    var url = window.location.href;
    var date = moment(ev.date);
    // find new date format
    if (ev.type == 'moveDate') {
        var newformat = 'day/' + date.format('YYYYMMDD') + '/';
        var dispdate = date.format('MMMM Do YYYY');
    }
    else if (ev.type == 'moveMonth') {
        var newformat = 'month/' + date.format('YYYYMM') + '/';
        var dispdate = date.format('MMMM YYYY');
    }
    else if (ev.type == 'moveYear') {
        var newformat = 'year/' + date.format('YYYY') + '/';
        var dispdate = date.format('YYYY');
    }
    // work through starting formats and proceed
    if (re_dayurl.test(url)) {
        var newurl = url.replace(re_dayurl, newformat);
    }
    else if (re_monthurl.test(url)) {
        var newurl = url.replace(re_dayurl, newformat);
    }
    else if (re_yearurl.test(url)) {
        var newurl = url.replace(re_dayurl, newformat);
    }
    else if (window.location.href ==
             document.getElementsByTagName('base')[0].href) {
        var newurl = window.location.href + newformat;
    }
    else {
        alert("ERROR: Cannot format new date. If the problem persists, please report this at https://github.com/gwpy/gwsumm/issues/");
    }
    /*$.ajax({
        type: 'HEAD',
        url : newurl,
        success: function(){window.location = newurl},
        error: function(){alert('No page found for ' + dispdate +
                                '. Please report unexpected problems to '+
                                'detchar+code@ligo.org.')}});
    */
    window.location = newurl;
}

// fix width of fixed navbar
function reset_width_on_resize() {
    $('#nav').width($("#nav").width());
}

// get basename of URL
function baseName(str) {
   base = new String(str).substring(str.lastIndexOf('/') + 1);
   return base; 
}

// add startsWith method
if (typeof String.prototype.startsWith != 'function') {
  String.prototype.startsWith = function (str){
    return this.slice(0, str.length) == str;
  };
}

function refreshImage() {
  var suffix = Math.floor((Math.random() * 10000) + 1);
  $('img').each(function() {
    var src = $(this).attr('src').split('?')[0];
    $(this).attr('src', src + "?" + suffix);
  });
}

function resizeFancyboxIframe() {
  var width = Math.min(1200, $(".fancybox-skin").width());
  // XXX: FIXME: this isn't finished yet
  if (width > document.body.clientWidth ) {
    $(".fancybox-iframe").width(width - 40);
  } else {
    $(".fancybox-iframe").width(width);
  }
  $(".fancybox-wrap").width(width + 30);

  // set heights as half width
  $(".fancybox-iframe").height(parseInt($(".fancybox-iframe").width() * 0.5));
  $(".fancybox-wrap").height(parseInt($(".fancybox-wrap").width() * 0.5));
}

/* ------------------------------------------------------------------------- */
/* Document ready and loaded                                                  */

// When document is ready, run this stuff:
$(window).load(function() {
    // define inter-IFO links
    var thisbase = document.getElementsByTagName('base')[0].href;
    $('a.ifo-switch').each(function() {
        var newbase = $(this).data('new-base') + '/';
        $(this).attr('href', window.location.href.replace(thisbase, newbase));
    });

    // define the calendar
    $('#calendar').datepicker({
        weekStart: 1,
        endDate: moment().utc().format('DD/MM/YYYY'),
        todayHighlight: true,
        todayBtn: "linked"
        }).on('moveDate', move_to_date)
              .on('moveMonth', move_to_date)
              .on('moveYear', move_to_date);

    // load correct run type
    if (location.hash.length > 1) {
        hash = location.hash.substring(1);
        path = location.pathname + hash + '.html';
        $('#state_' + hash).load_state(path);
    }

    // load the fancybox
    $(".fancybox").fancybox({
        nextEffect: 'none',
        prevEffect: 'none',
        width: 1200,
        iframe: {scrolling: 'no'},
        scrolling: 'no',
        beforeShow: function() {resizeFancyboxIframe()},
        helpers: {overlay: {locked: false}}
    });

    $('.dropdown-toggle').on('click', function() {
        var target = $(this).nextAll('.dropdown-menu');
        var dropleft = $(this).offset().left;
        var dropwidth = target.width();
        var left = $(window).width();
        var offright = (dropleft + dropwidth > left);
        if (offright) {
            target.addClass('pull-right');
        }
    });

    $(".fancybox-stamp").fancybox({
                    autoSize: false,
                    autoDimensions: false,
                    width: 1000,
                    height: 500,
                    fitToView: false,
                    showNavArrows: false,
                    padding: 0,
                    title: this.title,
                    href: $(this).attr('href'),
                    scrolling: 'no',
                    type: 'iframe'
                });



})

// When document is resized, fix the width of the sticky navbar
var resizeclock;
$(window).resize(function() {
    clearTimeout(resizeclock);
    $('#nav').width($('body').width());
    resizeclock = setTimeout(reset_width_on_resize, 100);
});
