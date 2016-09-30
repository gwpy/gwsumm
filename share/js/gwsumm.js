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

"use strict";

/* ------------------------------------------------------------------------- */
/* Calendar links                                                            */
var re_dayurl = new RegExp('day\/\\d{8}\/');
var re_monthurl = new RegExp('month\/\\d{6}\/');
var re_yearurl = new RegExp('year\/\\d{4}\/');

function findDateFormat() {
  var url = window.location.href;
  if (re_dayurl.test(url)) {
    return 'day';
  }
  if (re_monthurl.test(url)) {
    return 'month';
  }
  if (re_yearurl.test(url)) {
    return 'year';
  }
  return undefined;
}

function getPageDate() {
  var dateformat = findDateFormat();
  if (dateformat == 'day') {
    var datestring =
      String(re_dayurl.exec(window.location.href)).split('/')[1];
    return moment(datestring, 'YYYYMMDD');
  }
  if (dateformat == 'month') {
    var datestring =
      String(re_monthurl.exec(window.location.href)).split('/')[1];
    return moment(datestring, 'YYYYMM');
  }
  if (dateformat == 'year') {
    var datestring =
      String(re_yearurl.exec(window.location.href)).split('/')[1];
    return moment(datestring, 'YYYY');
  }
  throw "Cannot parse date format '" + dateformat + "'";
}

function stepDate(step) {
  var dateformat = findDateFormat();
  if (!dateformat) {
    return;
  }
  var date = getPageDate();
  var newdate = date.add(dateformat, step);
  if (dateformat == 'day') {
      move_to_date({type: 'changeDate', date: newdate});
  } else if (dateformat == 'month') {
      move_to_date({type: 'changeMonth', date: newdate});
  } else if (dateformat == 'year') {
      move_to_date({type: 'changeYear', date: newdate});
  }
}

// Load a given 'state'
$.fn.load_state = function loadState(page) {
  if ($(this).attr('id') == undefined) {
    return;
  }
  $('#main').load(page);
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
  if (ev.type == 'changeDate') {
    var newformat = 'day/' + date.format('YYYYMMDD') + '/';
    var dispdate = date.format('MMMM Do YYYY');
  }
  else if (ev.type == 'changeMonth') {
    var newformat = 'month/' + date.format('YYYYMM') + '/';
    var dispdate = date.format('MMMM YYYY');
  }
  else if (ev.type == 'changeYear') {
    var newformat = 'year/' + date.format('YYYY') + '/';
    var dispdate = date.format('YYYY');
  }
  // work through starting formats and proceed
  if (re_dayurl.test(url)) {
    var newurl = url.replace(re_dayurl, newformat);
  } else if (re_monthurl.test(url)) {
    var newurl = url.replace(re_dayurl, newformat);
  } else if (re_yearurl.test(url)) {
    var newurl = url.replace(re_dayurl, newformat);
  } else if (window.location.href ==
               document.getElementsByTagName('base')[0].href) {
    var newurl = window.location.href + newformat;
  } else {
    alert("ERROR: Cannot format new date. If the problem persists, please report this at https://github.com/gwpy/gwsumm/issues/");
  }
  window.location = newurl;
}

// fix width of fixed navbar
function reset_width_on_resize() {
  $('#nav').width($("#nav").width());
}

// resize fancybox iframe to 'normal' proportions
function resizeFancyboxIframe() {
  var width = Math.min(1200, $(".fancybox-skin").width());
  var height = (width - 40) * 0.5;
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

// shorten date in calendar if very small screen
function shortenDate() {
  var $calendar = $('#calendar');
  var date_ = moment($calendar.data('date'),
                     $calendar.data('date-format').toUpperCase());
  if ($calendar.html().startsWith('Calendar')) {  // don't break non-dates
    return;
  }
  if ($(document).width() < 400 ) {  // print shortened month name
    $calendar.html(date_.format('MMM D YYYY'));
  } else {  // print full month name
    $calendar.html(' ' + date_.format('MMMM D YYYY') + ' <b class="caret"></b>');
  }
}

/* ------------------------------------------------------------------------- */
/* Document ready and loaded                                                  */

// When document is ready, run this stuff:
$(window).load(function() {

  // shorten the date
  shortenDate();

  // define inter-IFO links
  var thisbase = document.getElementsByTagName('base')[0].href;
  $('[data-new-base]').each(function() {
    var newbase = $(this).data('new-base') + '/';
    $(this).attr('href', window.location.href.replace(thisbase, newbase));
  });

  // define the calendar
  $('#calendar').datepicker({
    weekStart: 1,
    endDate: moment().utc().format('DD/MM/YYYY'),
    todayHighlight: true,
    todayBtn: "linked"
  }).on('changeDate', move_to_date);

  // load correct run type
  if (location.hash.length > 1) {
    var hash = location.hash.substring(1);
    var path = location.pathname + hash + '.html';
    $('#state_' + hash).load_state(path);
  }

  // load the fancybox
  $(".fancybox").fancybox({nextEffect: 'none',
                           prevEffect: 'none',
                           width: 1200,
                           iframe: {scrolling: 'no'},
                           scrolling: 'no',
                           beforeShow: function() {resizeFancyboxIframe()},
                           helpers: {overlay: {locked: false}}
  });

  // custom fancybox for stamp-pem bokeh plot
  $(".fancybox-stamp").fancybox({width: 1000,
                                 height: 500,
                                 showNavArrows: false,
                                 padding: 0,
                                 title: this.title,
                                 href: $(this).attr('href'),
                                 type: 'iframe'
  });

  // reposition dropdown if too scrolling off the screen
  $('.dropdown-toggle').on('click', function() {
    // if page width is small, no-operation
    if ($(document).width() < 992 ) {
      return;
    }
    // otherwise add pull-right
    var target = $(this).nextAll('.dropdown-menu');
    var dropleft = $(this).offset().left;
    var dropwidth = target.width();
    var left = $(window).width();
    var offright = (dropleft + dropwidth > left);
    if (offright) {
      target.addClass('pull-right');
    }
  });

});

$(window).resize(function() {
  // set short month date
  shortenDate();
});
